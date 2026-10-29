-- | Render a QCircuit as a quantikz LaTeX environment.
--
-- Generates circuit diagrams.
-- 1-qubit QDiagonalAt [phi0, phi1] is rendered as P(phi1 - phi0);
-- global phase (phi0) is silently dropped.
-- n-qubit QDiagonalAt (n > 1) is rendered as a named box gate.
-- QGPhase (global phase) is omitted from the diagram.
module Quantikz where

import Numeric.Natural (Natural)
import Data.List (findIndex, intercalate)
import Data.Maybe (listToMaybe, mapMaybe)
import Text.Printf (printf)

import Syntax (QCircuit(..), QGate(..), gateQubits)

--------------------------------------------------------------------------------
-- Circuit width

circuitWidth :: QCircuit -> Natural
circuitWidth (QCircuit gates) =
  case concatMap gateQubits gates of
    [] -> 1
    qs -> maximum qs + 1

--------------------------------------------------------------------------------
-- Diagonal expansion

-- | Expand 1-qubit QDiagonalAt to QP, drop global phase.
--   Multi-qubit diagonals are left as-is (rendered as a box).
expandDiagonals :: [QGate] -> [QGate]
expandDiagonals = concatMap expand
  where
    expand (QDiagonalAt off 1 [phi0, phi1]) = [QP (phi1 - phi0) off]
    expand (QGPhase _) = []
    expand g           = [g]

--------------------------------------------------------------------------------
-- Column layout (greedy earliest-fit)

layoutGates :: [QGate] -> [[QGate]]
layoutGates = foldl insertGate []
  where
    insertGate cols g =
      case findIndex (canFit g) cols of
        Just i  -> [ if j == i then col ++ [g] else col
                   | (j, col) <- zip [0..] cols ]
        Nothing -> cols ++ [[g]]

    canFit g col =
      null [ q | q <- gateQubits g
               , q `elem` concatMap gateQubits col ]

--------------------------------------------------------------------------------
-- Per-cell rendering

renderCell :: Natural -> QGate -> String
renderCell q g = case g of
  QH i   | i == q -> "\\gate{H}"
  QS i   | i == q -> "\\gate{S}"
  QX i   | i == q -> "\\gate{X}"
  QP th i
    | i == q -> "\\gate{P(" ++ formatAngle th ++ ")}"
  QCX c t
    | c == q -> "\\ctrl{"  ++ relOff c t ++ "}"
    | t == q -> "\\targ{}"
  QSWAP a b
    | a == q -> "\\swap{"  ++ relOff a b ++ "}"
    | b == q -> "\\targX{}"
  QCP th c t
    | c == q -> "\\ctrl{"  ++ relOff c t ++ "}"
    | t == q -> "\\gate{P(" ++ formatAngle th ++ ")}"
  QDiagonalAt off n _phis
    | off == q ->
        "\\gate[wires=" ++ show n ++ "]{\\mathrm{D}}"
    | q > off && q < off + fromIntegral n ->
        ""
  _ -> "\\qw"
  where
    relOff a b = show (fromIntegral b - fromIntegral a :: Int)

--------------------------------------------------------------------------------
-- Full rendering

renderQuantikz :: QCircuit -> String
renderQuantikz qc =
  "\\begin{quantikz}[column sep=small,row sep={24pt,between origins}]\n"
  ++ intercalate " \\\\\n" (map wireRow [0 .. nq - 1])
  ++ "\n\\end{quantikz}"
  where
    nq   = circuitWidth qc
    cols = layoutGates (expandDiagonals (getQCircuit qc))

    wireRow q =
      "\\lstick{$q_{" ++ show q ++ "}$}"
      ++ concatMap (\col -> " & " ++ cellForWire q col) cols
      ++ " & \\qw"

    cellForWire q col =
      case mapMaybe (guardGate q) col of
        c : _ -> c
        []    -> "\\qw"

    guardGate q g
      | q `elem` gateQubits g = Just (renderCell q g)
      | otherwise              = Nothing

--------------------------------------------------------------------------------
-- Angle formatting (tries nice pi-fractions, falls back to decimal)

formatAngle :: Double -> String
formatAngle th
  | abs th < 1e-10 = "0"
  | otherwise      = case detectPiFraction th of
      Just s  -> s
      Nothing -> printf "%.4f" th

detectPiFraction :: Double -> Maybe String
detectPiFraction th = listToMaybe
  [ fmtFraction n d
  | d <- [1, 2, 4, 8, 16, 32, 64, 128, 256] :: [Int]
  , let ratio = th / pi
  , let n     = round (ratio * fromIntegral d) :: Int
  , n /= 0
  , abs (fromIntegral n / fromIntegral d - ratio) < 1e-6
  ]

fmtFraction :: Int -> Int -> String
fmtFraction n d
  | d == 1, n ==  1 = "\\pi"
  | d == 1, n == -1 = "{-}\\pi"
  | d == 1           = show n ++ "\\pi"
  | n ==  1          = "\\frac{\\pi}{" ++ show d ++ "}"
  | n == -1          = "{-}\\frac{\\pi}{" ++ show d ++ "}"
  | n  >  0          = "\\frac{" ++ show n ++ "\\pi}{" ++ show d ++ "}"
  | otherwise        = "{-}\\frac{" ++ show (abs n) ++ "\\pi}{" ++ show d ++ "}"
