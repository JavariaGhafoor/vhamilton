-- | Discrete noiseless VQE on the two-qubit DMANH+ chemistry Hamiltonian.
--
-- The target is the full dyadic Hamiltonian H_k, including the identity term.
--
-- The ansatz is a three-parameter real-state preparation family:
--
--   G1 = YI
--   G2 = 1/2 (IY + ZY) = |0><0| \otimes Y
--   G3 = 1/2 (IY - ZY) = |1><1| \otimes Y
--
-- Each parameter alpha_j = q_j / 2^r is fully discrete.
--
-- By default the classical optimizer uses a dyadic multiresolution beam search:
--
--   1. Exhaustively search a small coarse grid at r0 = min(3,r).
--   2. Keep the best B=3 coarse candidates.
--   3. On refinement r -> r+1, map each retained q to 2q, which represents
--      exactly the same physical alpha.
--   4. Search the +/-1 neighborhood of each doubled beam point.
--   5. Keep the best B candidates and repeat until the requested r.
--
-- This stays completely discrete while avoiding the 2^(3(r-1)) target-grid
-- brute-force cost.  A tighter coordinate mode replaces each full 3^p local
-- neighborhood by repeated axial +/-1 moves (about 2p points per beam center
-- per round).
--
-- The analytic mode goes one step further.  For the present three generators,
-- the energy as a function of one parameter (with the other two fixed) is a
-- trigonometric polynomial of degree at most two.  A constant number of VQE
-- energy samples therefore determines that one-dimensional energy curve.  We
-- fit the curve, minimize it classically, and jump directly to the target
-- dyadic grid rather than visiting every intermediate precision r.  This is a
-- dyadic, finite-grid analogue of Rotosolve / analytic coordinate
-- minimization.
--
-- The analytic-direct mode removes the remaining coarse exhaustive seed.  It
-- deterministically initializes every parameter at the balanced interior angle
-- beta_j = pi/4, i.e. alpha_j = 1/8 and q_j = 2^(r-3) for r >= 3.  For this
-- ansatz that prepares the equal-amplitude state |++>, so both conditional
-- branches are populated from the start.  The initialization needs no search
-- and is O(p) in the number of parameters; the mode then applies the same
-- analytic coordinate minimization directly at target precision r.  The old
-- exhaustive search remains available as a validation mode.
--
-- Usage:
--   cabal run chem-vqe
--   cabal run chem-vqe -- 10 5 0
--   cabal run chem-vqe -- 10 5 0 refine
--   cabal run chem-vqe -- 10 5 0 coordinate
--   cabal run chem-vqe -- 10 5 0 analytic
--   cabal run chem-vqe -- 10 5 0 analytic-direct
--   cabal run chem-vqe -- 10 5 0 exhaustive
--
-- Positional arguments are: Hamiltonian precision k, parameter precision r,
-- Trotter refinement n, and optional search mode.  Defaults are
-- k=10, r=3, n=0, mode=refine.
module Main where

import Control.Monad (foldM)
import Data.List (maximumBy, minimumBy, nub, sortBy)
import Data.Ord (comparing)
import Numeric.Natural (Natural)
import System.Environment (getArgs)
import System.IO (hFlush, stdout)
import Text.Printf (printf)
import Text.Read (readMaybe)

import ExamplesChem
  ( ChemDecomp(..)
  , fetchDecomp
  , chemProcWithIdentity
  , pauliParamProc
  , pauliSignedParamProc
  , identityDyadic
  )
import Syntax (Proc, (∥))
import VQE
  ( VBlock(..)
  , VAnsatz(..)
  , VQEResult(..)
  , energyForParams
  , gridSearch
  )

--------------------------------------------------------------------------------
-- Ansatz generators expressed directly in VHamilton

-- Search precision r affects only the external valuation q/2^r.
-- Atomic δ_k^alpha_j and -δ_k^alpha_j carry all variational dependence.

-- | Rotation of the first/high qubit: alpha1 * YI.
gFirst :: Proc
gFirst = pauliParamProc "YI" 0 "alpha1"

-- | Conditional Y on the second/low qubit when the first qubit is |0>:
-- alpha2/2 * (IY + ZY).
gCond0 :: Proc
gCond0 =
  pauliParamProc "IY" 1 "alpha2"
  ∥ pauliParamProc "ZY" 1 "alpha2"

-- | Conditional Y on the second/low qubit when the first qubit is |1>:
-- alpha3/2 * (IY - ZY).
gCond1 :: Proc
gCond1 =
  pauliParamProc "IY" 1 "alpha3"
  ∥ pauliSignedParamProc "ZY" 1 False "alpha3"

--------------------------------------------------------------------------------
-- CLI helpers

data SearchMode = Refine | Coordinate | Analytic | AnalyticDirect | Exhaustive
  deriving (Eq, Show)

readNatural :: String -> String -> Natural
readNatural name s =
  case readMaybe s :: Maybe Integer of
    Just x | x >= 0 -> fromInteger x
    _ -> error (name ++ " must be a non-negative integer")

readMode :: String -> SearchMode
readMode s =
  case s of
    "refine"     -> Refine
    "coordinate" -> Coordinate
    "coord"      -> Coordinate
    "analytic"        -> Analytic
    "rotosolve"       -> Analytic
    "analytic-direct" -> AnalyticDirect
    "direct"          -> AnalyticDirect
    "balanced"        -> AnalyticDirect
    "exhaustive"      -> Exhaustive
    _ -> error "search mode must be 'refine', 'coordinate' (or 'coord'), 'analytic' (or 'rotosolve'), 'analytic-direct' (or 'direct'/'balanced'), or 'exhaustive'"

parseArgs :: [String] -> (Natural, Natural, Natural, SearchMode)
parseArgs [] = (10, 3, 0, Refine)
parseArgs [k] = (readNatural "k" k, 3, 0, Refine)
parseArgs [k, r] = (readNatural "k" k, readNatural "r" r, 0, Refine)
parseArgs [k, r, n] =
  (readNatural "k" k, readNatural "r" r, readNatural "n" n, Refine)
parseArgs [k, r, n, mode] =
  ( readNatural "k" k
  , readNatural "r" r
  , readNatural "n" n
  , readMode mode
  )
parseArgs _ = error "usage: chem-vqe [k [r [n [refine|coordinate|analytic|analytic-direct|exhaustive]]]]"

-- We intentionally search alpha in [0, 1/2) for each parameter.  At precision
-- r this is q = 0,...,2^(r-1)-1.
parameterGrid :: Natural -> [Integer]
parameterGrid r
  | r < 1 = error "parameter precision r must be at least 1"
  | otherwise = [0 .. upper]
  where
    upper = 2 ^ (fromIntegral (r - 1) :: Int) - 1

blocksAt :: Natural -> Natural -> VAnsatz
blocksAt r n = VAnsatz
  { vaParamNames = ["alpha1", "alpha2", "alpha3"]
  , vaPrecision = r
  , vaBlocks =
      [ VBlock gFirst n
      , VBlock gCond0 n
      , VBlock gCond1 n
      ]
  }

--------------------------------------------------------------------------------
-- Dyadic multiresolution beam search

-- The coarse r=3 grid has only 4^3 = 64 candidates.  Retaining three seeds is
-- enough to preserve the relevant branch ambiguity in this benchmark while
-- still keeping each later refinement tiny.
coarsePrecision :: Natural
coarsePrecision = 3

beamWidth :: Int
beamWidth = 3

neighborRadius :: Integer
neighborRadius = 1

-- | All valid +/-radius neighbors of q at precision r, including q itself.
neighborhood :: Natural -> Integer -> [Integer] -> [[Integer]]
neighborhood r radius qs =
  sequence (map axis qs)
  where
    upper = 2 ^ (fromIntegral (r - 1) :: Int) - 1
    axis q =
      [ x
      | x <- [q - radius .. q + radius]
      , x >= 0
      , x <= upper
      ]

-- | Evaluate a finite candidate list with a live one-line progress display.
-- Returns results in evaluation order and the best result.
evaluateCandidates
  :: String
  -> Proc
  -> VAnsatz
  -> [[Integer]]
  -> IO ([VQEResult], VQEResult)
evaluateCandidates label target blocks candidates = do
  if null candidates
    then error "evaluateCandidates: empty candidate list"
    else pure ()

  let total = length candidates
      step
        :: (Int, [VQEResult], Maybe VQEResult)
        -> [Integer]
        -> IO (Int, [VQEResult], Maybe VQEResult)
      step (i0, acc, maybeBest) qs = do
        result <- energyForParams target blocks qs
        let i = i0 + 1
            best =
              case maybeBest of
                Nothing -> result
                Just old
                  | vrEnergy result < vrEnergy old -> result
                  | otherwise                      -> old
            pct = 100.0 * fromIntegral i / fromIntegral total :: Double

        printf
          "\r%s: %d/%d (%5.1f%%)  current q=%s  E=%.12f  best q=%s  best E=%.12f"
          label
          i
          total
          pct
          (show (vrParams result))
          (vrEnergy result)
          (show (vrParams best))
          (vrEnergy best)
        hFlush stdout
        pure (i, result : acc, Just best)

  (_, reversedResults, maybeBest) <-
    foldM step (0, [], Nothing) candidates
  putStrLn ""

  case maybeBest of
    Nothing -> error "evaluateCandidates: no results"
    Just best -> pure (reverse reversedResults, best)

-- | Coarse-to-fine discrete beam refinement.
--
-- At each new precision level, q -> 2q preserves alpha=q/2^r exactly.  The
-- newly introduced nearest grid points are then 2q +/- 1 on each coordinate.
-- Keeping a small beam protects against choosing a single coarse branch that
-- later turns out not to contain the fine-grid optimum.
dyadicRefinementSearch
  :: Proc
  -> Natural                 -- target precision r
  -> Natural                 -- Trotter refinement n
  -> IO (VQEResult, Integer) -- best result and number of energy evaluations
dyadicRefinementSearch target targetR n
  | targetR < 1 = error "parameter precision r must be at least 1"
  | otherwise = do
      let startR = min coarsePrecision targetR
          initialGrid = parameterGrid startR
          initialCandidates = sequence (replicate 3 initialGrid)

      putStrLn $ "Coarse search at r=" ++ show startR
      (initialResults, initialBest) <-
        evaluateCandidates
          ("r=" ++ show startR)
          target
          (blocksAt startR n)
          initialCandidates

      let initialBeam = take beamWidth (sortBy (comparing vrEnergy) initialResults)
          initialCount = toInteger (length initialCandidates)

      printf
        "  retained beam: %s\n"
        (show (map vrParams initialBeam))
      printf
        "  coarse best:   q=%s  E=%.12f\n\n"
        (show (vrParams initialBest))
        (vrEnergy initialBest)

      if startR == targetR
        then pure (initialBest, initialCount)
        else refineLevels (startR + 1) initialBeam initialCount
  where
    refineLevels
      :: Natural
      -> [VQEResult]
      -> Integer
      -> IO (VQEResult, Integer)
    refineLevels r beam totalEvals = do
      let doubled = map (map (* 2) . vrParams) beam
          candidates = nub (concatMap (neighborhood r neighborRadius) doubled)

      putStrLn $ "Refining to r=" ++ show r
      putStrLn $ "  inherited centers = " ++ show doubled
      putStrLn $ "  unique neighbors   = " ++ show (length candidates)

      (results, best) <-
        evaluateCandidates
          ("r=" ++ show r)
          target
          (blocksAt r n)
          candidates

      let nextBeam = take beamWidth (sortBy (comparing vrEnergy) results)
          totalEvals' = totalEvals + toInteger (length candidates)

      printf
        "  retained beam: %s\n"
        (show (map vrParams nextBeam))
      printf
        "  level best:    q=%s  E=%.12f\n\n"
        (show (vrParams best))
        (vrEnergy best)

      if r == targetR
        then pure (best, totalEvals')
        else refineLevels (r + 1) nextBeam totalEvals'

--------------------------------------------------------------------------------
-- Coordinate-wise dyadic multiresolution beam search

-- | Maximum number of axial local-refinement rounds at one dyadic precision.
-- In this benchmark convergence usually occurs in only a few rounds; this is
-- only a safeguard against pathological landscapes.
coordinateMaxRounds :: Int
coordinateMaxRounds = 20

-- | Axial +/-1 neighbors only.  Unlike 'neighborhood', this does not form the
-- Cartesian product across coordinates, so there are at most 2p new points
-- for p parameters rather than 3^p points.
axialNeighbors :: Natural -> [Integer] -> [[Integer]]
axialNeighbors r qs =
  [ replaceAt j x qs
  | (j, q) <- zip [0 ..] qs
  , x <- [q - 1, q + 1]
  , x >= 0
  , x <= upper
  ]
  where
    upper = 2 ^ (fromIntegral (r - 1) :: Int) - 1

    replaceAt j x xs =
      take j xs ++ [x] ++ drop (j + 1) xs

-- | Re-label a result from precision r as the exactly identical physical
-- point at precision r+1.  q/2^r = (2q)/2^(r+1), so no quantum energy
-- evaluation is needed for the inherited center.
inheritAtFinerPrecision :: VQEResult -> VQEResult
inheritAtFinerPrecision result =
  result { vrParams = map (* 2) (vrParams result) }

-- | Look up a cached result by its integer parameter tuple.
lookupCached :: [Integer] -> [VQEResult] -> Maybe VQEResult
lookupCached _ [] = Nothing
lookupCached qs (x:xs)
  | vrParams x == qs = Just x
  | otherwise        = lookupCached qs xs

-- | Keep the best B cached points.  The cache itself contains unique q-tuples
-- because candidates are filtered before evaluation.
bestBeam :: [VQEResult] -> [VQEResult]
bestBeam = take beamWidth . sortBy (comparing vrEnergy)

-- | Refine one precision level using only coordinate-axis moves.  Starting
-- from the inherited beam centers, repeatedly evaluate unseen +/-1 axial
-- neighbors, add them to a cache, and update the beam.  Stop when the beam no
-- longer changes or after 'coordinateMaxRounds' rounds.
coordinateRefineLevel
  :: Proc
  -> Natural
  -> Natural
  -> [VQEResult]
  -> IO ([VQEResult], Integer, Int)
coordinateRefineLevel target r n inherited =
  go 1 inherited inherited 0
  where
    go
      :: Int
      -> [VQEResult]  -- current beam
      -> [VQEResult]  -- all results cached at this level
      -> Integer      -- number of new quantum energy evaluations
      -> IO ([VQEResult], Integer, Int)
    go roundNo beam cache evalCount
      | roundNo > coordinateMaxRounds = do
          putStrLn $ "  coordinate safeguard reached " ++ show coordinateMaxRounds ++ " rounds"
          pure (bestBeam cache, evalCount, coordinateMaxRounds)
      | otherwise = do
          let proposed = nub (concatMap (axialNeighbors r . vrParams) beam)
              unseen =
                [ qs
                | qs <- proposed
                , case lookupCached qs cache of
                    Nothing -> True
                    Just _  -> False
                ]

          if null unseen
            then do
              putStrLn $ "  coordinate round " ++ show roundNo ++ ": no unseen axial neighbors"
              pure (bestBeam cache, evalCount, roundNo - 1)
            else do
              putStrLn $ "  coordinate round " ++ show roundNo
              putStrLn $ "    beam centers       = " ++ show (map vrParams beam)
              putStrLn $ "    new axial points   = " ++ show (length unseen)

              (newResults, _) <-
                evaluateCandidates
                  ("r=" ++ show r ++ " coord " ++ show roundNo)
                  target
                  (blocksAt r n)
                  unseen

              let cache' = newResults ++ cache
                  nextBeam = bestBeam cache'
                  evalCount' = evalCount + toInteger (length unseen)
                  stable = map vrParams nextBeam == map vrParams beam

              printf
                "    retained beam      = %s\n"
                (show (map vrParams nextBeam))
              printf
                "    round best         = q=%s  E=%.12f\n"
                (show (vrParams (head nextBeam)))
                (vrEnergy (head nextBeam))

              if stable
                then do
                  putStrLn "    beam unchanged; coordinate level converged"
                  pure (nextBeam, evalCount', roundNo)
                else go (roundNo + 1) nextBeam cache' evalCount'

-- | Coarse-to-fine coordinate beam refinement.
--
-- Complexity in energy evaluations after the coarse seed is approximately
-- O(B * p * L * r), where B is beam width, p is parameter count, and L is the
-- number of local axial rounds per refinement level.  This replaces the 3^p
-- neighborhood factor used by the full-neighborhood 'refine' mode.
coordinateDyadicSearch
  :: Proc
  -> Natural                 -- target precision r
  -> Natural                 -- Trotter refinement n
  -> IO (VQEResult, Integer) -- best result and number of actual evaluations
coordinateDyadicSearch target targetR n
  | targetR < 1 = error "parameter precision r must be at least 1"
  | otherwise = do
      let startR = min coarsePrecision targetR
          initialGrid = parameterGrid startR
          initialCandidates = sequence (replicate 3 initialGrid)

      putStrLn $ "Coarse search at r=" ++ show startR
      (initialResults, initialBest) <-
        evaluateCandidates
          ("r=" ++ show startR)
          target
          (blocksAt startR n)
          initialCandidates

      let initialBeam = bestBeam initialResults
          initialCount = toInteger (length initialCandidates)

      printf
        "  retained beam: %s\n"
        (show (map vrParams initialBeam))
      printf
        "  coarse best:   q=%s  E=%.12f\n\n"
        (show (vrParams initialBest))
        (vrEnergy initialBest)

      if startR == targetR
        then pure (initialBest, initialCount)
        else refineLevels (startR + 1) initialBeam initialCount
  where
    refineLevels
      :: Natural
      -> [VQEResult]
      -> Integer
      -> IO (VQEResult, Integer)
    refineLevels r previousBeam totalEvals = do
      let inherited = map inheritAtFinerPrecision previousBeam

      putStrLn $ "Coordinate refining to r=" ++ show r
      putStrLn $ "  inherited centers = " ++ show (map vrParams inherited)
      putStrLn $ "  inherited centers reuse their previous energies (0 new evaluations)"

      (nextBeam, newEvals, rounds) <-
        coordinateRefineLevel target r n inherited

      let best = head nextBeam
          totalEvals' = totalEvals + newEvals

      printf "  coordinate rounds  = %d\n" rounds
      printf "  new evaluations    = %d\n" newEvals
      printf
        "  level best:         q=%s  E=%.12f\n\n"
        (show (vrParams best))
        (vrEnergy best)

      if r == targetR
        then pure (best, totalEvals')
        else refineLevels (r + 1) nextBeam totalEvals'

--------------------------------------------------------------------------------
-- Analytic dyadic coordinate minimization

-- For the current ansatz, fixing two parameters and varying the third gives
-- a low-degree trigonometric polynomial in
--
--   phi = 2*pi*alpha = pi*q / 2^(r-1).
--
-- G1 = YI has spectrum {-1,+1}, so only the frequency-2 harmonic can appear:
--
--   E_1(phi) = c + a2 cos(2 phi) + b2 sin(2 phi).
--
-- G2 and G3 have spectrum {-1,0,+1}, so frequency differences 1 and 2 can
-- appear:
--
--   E_j(phi) = c + a1 cos(phi) + b1 sin(phi)
--                + a2 cos(2 phi) + b2 sin(2 phi).
--
-- Hence three energy samples determine the G1 coordinate curve and five
-- samples determine the G2/G3 coordinate curves.  Once the curve is fitted,
-- its minimum is found with only classical arithmetic and then projected back
-- to the target dyadic grid.  The expensive VQE-energy query count is thus
-- independent of the requested precision r, apart from the fixed coarse seed.
-- AnalyticDirect below removes that seed search as well.

data TrigModel = TrigModel
  { tmC  :: Double
  , tmA1 :: Double
  , tmB1 :: Double
  , tmA2 :: Double
  , tmB2 :: Double
  } deriving (Eq, Show)

analyticMaxSweeps :: Int
analyticMaxSweeps = 12

analyticScanBins :: Int
analyticScanBins = 256

analyticGoldenIterations :: Int
analyticGoldenIterations = 70

analyticResidualTolerance :: Double
analyticResidualTolerance = 1.0e-7

gridSizeAt :: Natural -> Integer
gridSizeAt r
  | r < 1 = error "gridSizeAt: r must be at least 1"
  | otherwise = 2 ^ (fromIntegral (r - 1) :: Int)

replaceCoordinate :: Int -> Integer -> [Integer] -> [Integer]
replaceCoordinate j x xs =
  take j xs ++ [x] ++ drop (j + 1) xs

phaseForQ :: Integer -> Integer -> Double
phaseForQ m q =
  pi * fromIntegral q / fromIntegral m

trigValue :: TrigModel -> Double -> Double
trigValue model phi =
  tmC model
    + tmA1 model * cos phi
    + tmB1 model * sin phi
    + tmA2 model * cos (2 * phi)
    + tmB2 model * sin (2 * phi)

-- | Small Gauss-Jordan solver with partial pivoting.  It is used only for the
-- fixed 5x5 trigonometric fit, so keeping it local avoids another dependency.
solveLinear :: [[Double]] -> [Double] -> [Double]
solveLinear a b
  | null a = error "solveLinear: empty matrix"
  | length a /= length b = error "solveLinear: matrix/vector size mismatch"
  | any ((/= n) . length) a = error "solveLinear: matrix must be square"
  | otherwise = map last (foldl eliminateColumn augmented [0 .. n - 1])
  where
    n = length a
    augmented = zipWith (++) a (map (:[]) b)

    eliminateColumn :: [[Double]] -> Int -> [[Double]]
    eliminateColumn rows k =
      let candidates =
            [ (i, abs ((rows !! i) !! k))
            | i <- [k .. n - 1]
            ]
          (pivotIndex, pivotAbs) = maximumBy (comparing snd) candidates
      in if pivotAbs < 1.0e-12
           then error "solveLinear: singular trigonometric fit"
           else
             let rowsSwapped = swapRows k pivotIndex rows
                 pivotRow0 = rowsSwapped !! k
                 pivotValue = pivotRow0 !! k
                 pivotRow = map (/ pivotValue) pivotRow0
                 reduceRow i row
                   | i == k = pivotRow
                   | otherwise =
                       let factor = row !! k
                       in zipWith (-) row (map (* factor) pivotRow)
             in zipWith reduceRow [0 ..] rowsSwapped

    swapRows :: Int -> Int -> [[Double]] -> [[Double]]
    swapRows i j rows
      | i == j = rows
      | otherwise =
          [ choose k row
          | (k, row) <- zip [0 ..] rows
          ]
      where
        ri = rows !! i
        rj = rows !! j
        choose k row
          | k == i = rj
          | k == j = ri
          | otherwise = row

-- | Fit the exact one-coordinate trigonometric model for the present ansatz.
fitTrigModel
  :: Int       -- coordinate 0,1,2
  -> Integer   -- M = 2^(r-1)
  -> [Integer] -- sampled q_j values
  -> [Double]  -- sampled energies
  -> TrigModel
fitTrigModel j m sampleQs energies
  | j == 0 =
      case energies of
        [e0, eQuarter, eHalf] ->
          let c = (e0 + eHalf) / 2
              a2 = (e0 - eHalf) / 2
              b2 = eQuarter - c
          in TrigModel c 0 0 a2 b2
        _ -> error "fitTrigModel: G1 requires exactly three samples"
  | j == 1 || j == 2 =
      let row q =
            let phi = phaseForQ m q
            in [ 1
               , cos phi
               , sin phi
               , cos (2 * phi)
               , sin (2 * phi)
               ]
          coeffs = solveLinear (map row sampleQs) energies
      in case coeffs of
           [c, a1, b1, a2, b2] -> TrigModel c a1 b1 a2 b2
           _ -> error "fitTrigModel: internal coefficient-size error"
  | otherwise = error "fitTrigModel: coordinate index out of range"

-- | Sample locations used to reconstruct one coordinate curve.  All samples
-- remain inside the same allowed half-period q-domain [0,M).
analyticSampleQs :: Int -> Integer -> [Integer]
analyticSampleQs j m
  | j == 0 = [0, m `div` 4, m `div` 2]
  | j == 1 || j == 2 =
      [0, m `div` 4, m `div` 2, 3 * m `div` 4, 7 * m `div` 8]
  | otherwise = error "analyticSampleQs: coordinate index out of range"

-- | Golden-section minimization on a short bracket.  The surrounding coarse
-- scan identifies every visible local-minimum basin of this degree-2 model;
-- golden section then refines each basin without any VQE calls.
goldenMin :: TrigModel -> Double -> Double -> Double
goldenMin model a0 b0 = go analyticGoldenIterations a0 b0 c0 d0 fc0 fd0
  where
    rho = (sqrt 5 - 1) / 2
    c0 = b0 - rho * (b0 - a0)
    d0 = a0 + rho * (b0 - a0)
    fc0 = trigValue model c0
    fd0 = trigValue model d0

    go 0 a b _ _ _ _ = (a + b) / 2
    go it a b c d fc fd
      | fc < fd =
          let b' = d
              d' = c
              fd' = fc
              c' = b' - rho * (b' - a)
              fc' = trigValue model c'
          in go (it - 1) a b' c' d' fc' fd'
      | otherwise =
          let a' = c
              c' = d
              fc' = fd
              d' = a' + rho * (b - a')
              fd' = trigValue model d'
          in go (it - 1) a' b c' d' fc' fd'

-- | Candidate continuous minima on phi in [0,pi].  The upper endpoint is
-- included only as a boundary probe; q=M itself is not in the allowed grid,
-- and 'bestDyadicQ' maps it back to the nearest valid boundary point.
modelMinimumPhases :: TrigModel -> [Double]
modelMinimumPhases model =
  0 : pi : map refine localBrackets
  where
    nBins = analyticScanBins
    phi i = pi * fromIntegral i / fromIntegral nBins
    value i = trigValue model (phi i)
    localBrackets =
      [ (phi (i - 1), phi (i + 1))
      | i <- [1 .. nBins - 1]
      , value i <= value (i - 1)
      , value i <= value (i + 1)
      ]
    refine (a, b) = goldenMin model a b

-- | Project all continuous local minima (and both boundaries) to a small set
-- of nearby valid dyadic q values, then choose the lowest point according to
-- the fitted model.  The number of classical candidates is bounded by a small
-- constant because the model degree is fixed.
bestDyadicQ :: Integer -> TrigModel -> Integer
bestDyadicQ m model =
  minimumBy (comparing modelEnergy) candidates
  where
    upper = m - 1
    clamp q = max 0 (min upper q)
    around phi =
      let qContinuous = fromIntegral m * phi / pi
          qFloor = floor qContinuous :: Integer
          qCeil = ceiling qContinuous :: Integer
          qRound = round qContinuous :: Integer
      in map clamp [qFloor, qCeil, qRound - 1, qRound, qRound + 1]
    candidates = nub ([0, upper] ++ concatMap around (modelMinimumPhases model))
    modelEnergy q = trigValue model (phaseForQ m q)

-- | Evaluate candidate tuples while reusing any values already present in a
-- cache.  This is important in analytic mode because the same sample points
-- recur across coordinate fits and beam seeds.
evaluateWithCache
  :: String
  -> Proc
  -> VAnsatz
  -> [VQEResult]
  -> [[Integer]]
  -> IO ([VQEResult], [VQEResult], Integer)
evaluateWithCache label target blocks cache requested = do
  let uniqueRequested = nub requested
      unseen =
        [ qs
        | qs <- uniqueRequested
        , case lookupCached qs cache of
            Nothing -> True
            Just _  -> False
        ]

  newResults <-
    if null unseen
      then pure []
      else fst <$> evaluateCandidates label target blocks unseen

  let cache' = newResults ++ cache
      fetch qs =
        case lookupCached qs cache' of
          Just result -> result
          Nothing -> error "evaluateWithCache: requested point missing from cache"

  pure (map fetch requested, cache', toInteger (length unseen))

-- | Remove duplicate q-tuples while keeping the first result for each tuple.
uniqueResults :: [VQEResult] -> [VQEResult]
uniqueResults [] = []
uniqueResults (x:xs) =
  x : uniqueResults [y | y <- xs, vrParams y /= vrParams x]

-- | Perform one sequential analytic coordinate pass for one beam seed.
analyticPass
  :: Proc
  -> Natural
  -> Natural
  -> Int
  -> Int
  -> VQEResult
  -> [VQEResult]
  -> Integer
  -> IO (VQEResult, [VQEResult], Integer)
analyticPass target r n sweepNo seedNo startResult cache0 evalCount0 =
  go 0 startResult cache0 evalCount0
  where
    m = gridSizeAt r
    blocks = blocksAt r n

    go
      :: Int
      -> VQEResult
      -> [VQEResult]
      -> Integer
      -> IO (VQEResult, [VQEResult], Integer)
    go j current cache evalCount
      | j >= 3 = pure (current, cache, evalCount)
      | otherwise = do
          let qs = vrParams current
              sampleQs = analyticSampleQs j m
              sampleTuples = [replaceCoordinate j q qs | q <- sampleQs]
              fitLabel =
                "analytic sweep " ++ show sweepNo
                  ++ " seed " ++ show seedNo
                  ++ " q" ++ show (j + 1) ++ " fit"

          (sampleResults, cache1, addedFit) <-
            evaluateWithCache fitLabel target blocks cache sampleTuples

          let model = fitTrigModel j m sampleQs (map vrEnergy sampleResults)
              qBest = bestDyadicQ m model
              predictedTuple = replaceCoordinate j qBest qs
              jumpLabel =
                "analytic sweep " ++ show sweepNo
                  ++ " seed " ++ show seedNo
                  ++ " q" ++ show (j + 1) ++ " jump"

          (jumpResults, cache2, addedJump) <-
            evaluateWithCache jumpLabel target blocks cache1 [predictedTuple]

          let next = head jumpResults
              predictedEnergy = trigValue model (phaseForQ m qBest)
              residual = vrEnergy next - predictedEnergy
              evalCount' = evalCount + addedFit + addedJump

          printf
            "    q%d analytic jump: %d -> %d  E=%.12f  model residual=%+.3e\n"
            (j + 1)
            (qs !! j)
            qBest
            (vrEnergy next)
            residual

          if abs residual > analyticResidualTolerance
            then putStrLn
              "      WARNING: fitted harmonic model residual is larger than expected; analytic mode may not be valid for this generator/compiler setting."
            else pure ()

          go (j + 1) next cache2 evalCount'

-- | Direct-to-target analytic coordinate optimization.  A small r=3
-- exhaustive search supplies a branch-safe beam, whose physical parameters
-- are lifted exactly to the requested r in one step.  Each coordinate is then
-- minimized from a constant-size trigonometric fit.  Consequently the number
-- of expensive VQE energy evaluations is approximately
--
--   O(2^(p*r0) + B*p*L),
--
-- with fixed r0=3, beam width B, parameter count p=3, and coordinate-sweep
-- count L.  In particular, there is no multiplicative target-r factor in the
-- energy-query complexity.
analyticDyadicSearch
  :: Proc
  -> Natural
  -> Natural
  -> IO (VQEResult, Integer)
analyticDyadicSearch target targetR n
  | targetR < 1 = error "parameter precision r must be at least 1"
  | targetR <= coarsePrecision = do
      let grid = parameterGrid targetR
          candidates = sequence (replicate 3 grid)
      putStrLn $ "Analytic mode: target r=" ++ show targetR
        ++ " is no finer than r0=" ++ show coarsePrecision
        ++ "; using the exact small grid instead."
      (results, best) <-
        evaluateCandidates
          ("r=" ++ show targetR)
          target
          (blocksAt targetR n)
          candidates
      pure (best, toInteger (length results))
  | otherwise = do
      let startR = coarsePrecision
          initialGrid = parameterGrid startR
          initialCandidates = sequence (replicate 3 initialGrid)

      putStrLn $ "Coarse search at r=" ++ show startR
      (initialResults, initialBest) <-
        evaluateCandidates
          ("r=" ++ show startR)
          target
          (blocksAt startR n)
          initialCandidates

      let initialBeam = bestBeam initialResults
          liftFactor = 2 ^ (fromIntegral (targetR - startR) :: Int)
          liftedBeam =
            [ result { vrParams = map (* liftFactor) (vrParams result) }
            | result <- initialBeam
            ]
          initialCount = toInteger (length initialCandidates)

      printf "  retained coarse beam: %s\n" (show (map vrParams initialBeam))
      printf
        "  coarse best:          q=%s  E=%.12f\n"
        (show (vrParams initialBest))
        (vrEnergy initialBest)
      putStrLn $ "  direct lift to r=" ++ show targetR
      putStrLn $ "  lifted beam:          " ++ show (map vrParams liftedBeam)
      putStrLn "  lifted points reuse their coarse energies (0 new evaluations)"
      putStrLn ""

      sweep 1 liftedBeam liftedBeam initialCount
  where
    sweep
      :: Int
      -> [VQEResult]
      -> [VQEResult]
      -> Integer
      -> IO (VQEResult, Integer)
    sweep sweepNo beam cache totalEvals
      | sweepNo > analyticMaxSweeps = do
          putStrLn $ "Analytic safeguard reached " ++ show analyticMaxSweeps ++ " sweeps"
          pure (head (bestBeam cache), totalEvals)
      | otherwise = do
          putStrLn $ "Analytic coordinate sweep " ++ show sweepNo
          putStrLn $ "  beam starts = " ++ show (map vrParams beam)

          (updatedRev, cache', totalEvals') <-
            foldM
              (\(acc, cacheAcc, evalAcc) (seedNo, seed) -> do
                  putStrLn $ "  seed " ++ show seedNo ++ " start q=" ++ show (vrParams seed)
                  (updated, cacheNext, evalNext) <-
                    analyticPass
                      target
                      targetR
                      n
                      sweepNo
                      seedNo
                      seed
                      cacheAcc
                      evalAcc
                  putStrLn $ "  seed " ++ show seedNo ++ " end   q=" ++ show (vrParams updated)
                    ++ "  E=" ++ printf "%.12f" (vrEnergy updated)
                  pure (updated : acc, cacheNext, evalNext)
              )
              ([], cache, totalEvals)
              (zip [1 ..] beam)

          let updated = reverse updatedRev
              nextBeam = bestBeam (uniqueResults (updated ++ beam))
              stable = map vrParams nextBeam == map vrParams beam
              sweepNewEvals = totalEvals' - totalEvals

          putStrLn $ "  retained beam = " ++ show (map vrParams nextBeam)
          putStrLn $ "  sweep evaluations = " ++ show sweepNewEvals
          putStrLn $ "  sweep best = q=" ++ show (vrParams (head nextBeam))
            ++ "  E=" ++ printf "%.12f" (vrEnergy (head nextBeam))

          if stable
            then do
              putStrLn "  beam unchanged; analytic coordinate search converged"
              putStrLn ""
              pure (head nextBeam, totalEvals')
            else do
              putStrLn ""
              sweep (sweepNo + 1) nextBeam cache' totalEvals'

--------------------------------------------------------------------------------
-- Deterministic balanced-start analytic coordinate minimization

-- | Balanced interior dyadic parameter.  For r >= 3,
--
--     q / 2^r = 1/8
--
-- exactly, so beta = 2*pi*alpha = pi/4.  Applying this value to all three
-- generators prepares [1,1,1,1]/2 = |++> for the present ansatz.  In
-- particular, neither conditional branch starts with zero population.
balancedSeedQ :: Natural -> Integer
balancedSeedQ r
  | r < 3 = error "balancedSeedQ: exact pi/4 dyadic seed requires r >= 3"
  | otherwise = 2 ^ (fromIntegral (r - 3) :: Int)

balancedSeed :: Natural -> [Integer]
balancedSeed r = replicate 3 (balancedSeedQ r)

-- | Analytic coordinate optimization started directly from the deterministic
-- balanced interior point at the target precision.  There is no coarse
-- exhaustive search and no traversal through intermediate r values.
--
-- The start itself is constructed with O(p) integer work and zero VQE queries;
-- one energy query records its initial objective value.  Thereafter each
-- coordinate uses a constant-size exact trigonometric fit, so the expensive
-- objective-query complexity is O(p*L) for fixed harmonic degree, where L is
-- the number of analytic coordinate sweeps.  This is a heuristic global
-- optimizer: the polynomial initialization is guaranteed, but convergence to
-- the global VQE minimum is not guaranteed for a general nonconvex landscape.
analyticDirectSearch
  :: Proc
  -> Natural
  -> Natural
  -> IO (VQEResult, Integer)
analyticDirectSearch target targetR n
  | targetR < 1 = error "parameter precision r must be at least 1"
  | targetR < 3 = do
      let grid = parameterGrid targetR
          candidates = sequence (replicate 3 grid)
      putStrLn $ "Analytic-direct mode: r=" ++ show targetR
        ++ " cannot represent the balanced alpha=1/8 seed exactly;"
      putStrLn "using the exact tiny grid instead."
      (results, best) <-
        evaluateCandidates
          ("r=" ++ show targetR)
          target
          (blocksAt targetR n)
          candidates
      pure (best, toInteger (length results))
  | otherwise = do
      let q0 = balancedSeed targetR
          blocks = blocksAt targetR n

      putStrLn $ "Deterministic balanced start at target r=" ++ show targetR
      putStrLn "  beta_j  = pi/4"
      putStrLn "  alpha_j = 1/8"
      putStrLn $ "  q_j     = 2^(r-3) = " ++ show (balancedSeedQ targetR)
      putStrLn $ "  start q = " ++ show q0
      putStrLn "  prepared reference state = |++> = [1,1,1,1]/2"
      putStrLn "  initialization search cost = 0 candidates"

      (seedResults, seedBest) <-
        evaluateCandidates "balanced seed" target blocks [q0]

      putStrLn $ "  start energy = " ++ printf "%.12f" (vrEnergy seedBest)
      putStrLn ""

      sweepDirect 1 seedBest seedResults 1
  where
    sweepDirect
      :: Int
      -> VQEResult
      -> [VQEResult]
      -> Integer
      -> IO (VQEResult, Integer)
    sweepDirect sweepNo current cache totalEvals
      | sweepNo > analyticMaxSweeps = do
          putStrLn $ "Analytic-direct safeguard reached " ++ show analyticMaxSweeps ++ " sweeps"
          pure (current, totalEvals)
      | otherwise = do
          putStrLn $ "Analytic-direct coordinate sweep " ++ show sweepNo
          putStrLn $ "  start q = " ++ show (vrParams current)

          (updated, cache', totalEvals') <-
            analyticPass
              target
              targetR
              n
              sweepNo
              1
              current
              cache
              totalEvals

          let changed = vrParams updated /= vrParams current
              improved = vrEnergy updated <= vrEnergy current + 1.0e-10
              next = if improved then updated else current
              accepted = vrParams next /= vrParams current
              sweepNewEvals = totalEvals' - totalEvals

          putStrLn $ "  end q   = " ++ show (vrParams updated)
            ++ "  E=" ++ printf "%.12f" (vrEnergy updated)
          putStrLn $ "  sweep evaluations = " ++ show sweepNewEvals

          if not improved
            then do
              putStrLn "  fitted coordinate pass increased energy numerically; retaining previous point"
              putStrLn "  analytic-direct search converged"
              putStrLn ""
              pure (current, totalEvals')
            else if not changed || not accepted
              then do
                putStrLn "  q unchanged; analytic-direct coordinate search converged"
                putStrLn ""
                pure (next, totalEvals')
              else do
                putStrLn ""
                sweepDirect (sweepNo + 1) next cache' totalEvals'

--------------------------------------------------------------------------------

main :: IO ()
main = do
  args <- getArgs
  let (k, r, n, mode) = parseArgs args
      grid = parameterGrid r
      exhaustiveCount = toInteger (length grid) ^ (3 :: Int)

  d <- fetchDecomp "upper" k

  let target = chemProcWithIdentity d

  putStrLn "Discrete chemistry VQE: DMANH+ upper 4x4 block"
  putStrLn ""
  putStrLn $ "Hamiltonian dyadic precision k = " ++ show k
  putStrLn $ "parameter dyadic precision r   = " ++ show r
  putStrLn $ "Trotter refinement n           = " ++ show n
  putStrLn $ "search mode                    = " ++
    case mode of
      Refine     -> "dyadic multiresolution beam refinement (full 3^p neighborhood)"
      Coordinate -> "coordinate-wise dyadic beam refinement (axial 2p neighborhood)"
      Analytic       -> "analytic dyadic coordinate minimization (trigonometric fit + coarse beam seed)"
      AnalyticDirect -> "analytic-direct minimization (deterministic balanced interior seed)"
      Exhaustive     -> "full exhaustive grid"
  putStrLn $ "target-grid candidates         = " ++ show exhaustiveCount
  case mode of
    Refine -> do
      putStrLn $ "coarse precision r0            = " ++ show (min coarsePrecision r)
      putStrLn $ "beam width                     = " ++ show beamWidth
      putStrLn $ "neighbor radius                = +/-" ++ show neighborRadius
    Coordinate -> do
      putStrLn $ "coarse precision r0            = " ++ show (min coarsePrecision r)
      putStrLn $ "beam width                     = " ++ show beamWidth
      putStrLn $ "axial moves per center         = at most " ++ show (2 * 3 :: Int)
      putStrLn $ "max coordinate rounds/level    = " ++ show coordinateMaxRounds
    Analytic -> do
      putStrLn $ "coarse precision r0            = " ++ show (min coarsePrecision r)
      putStrLn $ "beam width                     = " ++ show beamWidth
      putStrLn "G1 fit samples                 = 3 energy values"
      putStrLn "G2/G3 fit samples              = 5 energy values each"
      putStrLn $ "max analytic coordinate sweeps = " ++ show analyticMaxSweeps
      putStrLn "precision traversal            = direct r0 -> target r"
    AnalyticDirect -> do
      putStrLn "initialization                  = deterministic balanced interior point"
      putStrLn "balanced angle                  = beta_j = pi/4, alpha_j = 1/8"
      if r >= 3
        then putStrLn $ "balanced integer q_j            = 2^(r-3) = " ++ show (balancedSeedQ r)
        else putStrLn "balanced integer q_j            = unavailable exactly for r < 3; tiny exhaustive fallback"
      putStrLn "G1 fit samples                 = 3 energy values"
      putStrLn "G2/G3 fit samples              = 5 energy values each"
      putStrLn $ "max analytic coordinate sweeps = " ++ show analyticMaxSweeps
      putStrLn "precision traversal            = none; optimize directly at target r"
      putStrLn "coarse exhaustive seed         = none"
    Exhaustive ->
      putStrLn $ "q-domain per parameter         = " ++ show grid
  putStrLn ""
  putStrLn "Ansatz generators:"
  putStrLn "  G1 = YI"
  putStrLn "  G2 = (IY + ZY)/2 = |0><0| x Y"
  putStrLn "  G3 = (IY - ZY)/2 = |1><1| x Y"
  putStrLn "  alpha_j = q_j / 2^r"
  putStrLn ""

  (best, evaluatedCount) <-
    case mode of
      Refine -> dyadicRefinementSearch target r n
      Coordinate -> coordinateDyadicSearch target r n
      Analytic -> analyticDyadicSearch target r n
      AnalyticDirect -> analyticDirectSearch target r n
      Exhaustive -> do
        let grids = replicate 3 grid
        (b, _) <- gridSearch target (blocksAt r n) grids
        pure (b, exhaustiveCount)

  let ident = identityDyadic d
      scaledVQEFull = vrEnergy best
      scaledHkGround = vrGroundEnergy best
      nonidentityVQE = scaledVQEFull - ident
      nonidentityGround = scaledHkGround - ident
      physicalVQE = scaledVQEFull / cdScale d
      physicalHkGround = scaledHkGround / cdScale d
      variationalGapScaled = scaledVQEFull - scaledHkGround
      variationalGapMh = variationalGapScaled / cdScale d
      hamDiscSignedMh = physicalHkGround - cdPhysicalGroundMh d
      totalSignedMh = physicalVQE - cdPhysicalGroundMh d
      reduction =
        if evaluatedCount == 0
          then 0.0
          else fromIntegral exhaustiveCount / fromIntegral evaluatedCount :: Double

  putStrLn "Best discrete VQE result"
  putStrLn $ "  best q                     = " ++ show (vrParams best)
  putStrLn $ "  full VQE energy (scaled)   = " ++ printf "%.12f" scaledVQEFull
  putStrLn $ "  exact ground of dyadic H_k (scaled) = " ++ printf "%.12f" scaledHkGround
  putStrLn $ "  nonidentity VQE component  = " ++ printf "%.12f" nonidentityVQE
  putStrLn $ "  nonidentity ground comp.   = " ++ printf "%.12f" nonidentityGround
  putStrLn $ "  variational gap (scaled)   = " ++ printf "%.12e" variationalGapScaled
  putStrLn ""

  putStrLn "Search cost"
  putStrLn $ "  energy evaluations used    = " ++ show evaluatedCount
  putStrLn $ "  exhaustive target grid     = " ++ show exhaustiveCount
  putStrLn $ "  evaluation reduction       = " ++ printf "%.1fx" reduction
  putStrLn ""

  putStrLn "Identity/physical-energy reconstruction"
  putStrLn $ "  scale                      = " ++ printf "%.12f" (cdScale d)
  putStrLn $ "  exact identity coeff (mH) = " ++ printf "%.12f" (cdIdentityCoeffMh d)
  putStrLn $ "  scaled identity (exact)    = " ++ printf "%.12f" (cdScaledIdentity d)
  putStrLn $ "  scaled identity (dyadic)   = " ++ printf "%.12f" ident
  putStrLn ""

  putStrLn "Physical energies (milliHartree)"
  putStrLn $ "  VQE estimate               = " ++ printf "%.12f" physicalVQE
  putStrLn $ "  exact ground of dyadic H_k = " ++ printf "%.12f" physicalHkGround
  putStrLn $ "  exact ground of physical H = " ++ printf "%.12f" (cdPhysicalGroundMh d)
  putStrLn $ "  variational gap            = " ++ printf "%.12f" variationalGapMh
  putStrLn $ "  H discretization shift     = " ++ printf "%+.12f" hamDiscSignedMh
  putStrLn $ "  total signed VQE error     = " ++ printf "%+.12f" totalSignedMh
  putStrLn ""
  putStrLn "Notes:"
  putStrLn "  * The VQE target includes the dyadic identity term; ansatz circuits do not need it."
  putStrLn "  * Refine, coordinate, analytic, and analytic-direct modes all return parameters on the requested dyadic grid."
  putStrLn "  * Coordinate mode tests axial +/-1 moves first, reducing the local 3^p factor to roughly 2p per round."
  putStrLn "  * Analytic mode fits the exact low-degree one-coordinate trigonometric energy curves for this ansatz and jumps directly to target r."
  putStrLn "  * Analytic mode has energy-query complexity approximately O(2^(p*r0) + B*p*L)."
  putStrLn "  * Analytic-direct removes the coarse search by using q_j=2^(r-3), the deterministic beta_j=pi/4 balanced seed."
  putStrLn "  * For fixed harmonic degree its energy-query complexity is approximately O(p*L), independent of target precision r."
  putStrLn "  * The balanced-start construction is polynomial; global convergence of the nonconvex VQE search is still heuristic."
  putStrLn "  * All non-exhaustive search modes remain heuristic as global multi-parameter optimizers."
  putStrLn "  * Use 'exhaustive' mode at feasible r to validate the refinement result."
  putStrLn "  * G2 and G3 contain commuting Pauli terms, so n=0 has no Trotter error here."
  putStrLn "  * This run is exact-statevector/noiseless; there is no shot noise yet."
