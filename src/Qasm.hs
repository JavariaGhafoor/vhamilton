-- | Emit QCircuit as OpenQASM 3 (§Implementation, fig:circuits).
--
-- Final stage of the compilation pipeline:  QCircuit → OpenQASM 3 string.
-- Multi-qubit QDiagonalAt gates are lowered to elementary gates by
-- scripts.compiler.diag_to_qasm3, which uses Qiskit's DiagonalGate and transpiles.
module Qasm where

import Numeric.Natural (Natural)
import GHC.Generics (Generic)
import qualified Data.Aeson as A
import qualified Data.ByteString.Lazy.Char8 as BL
import System.Exit (ExitCode(..))
import System.Process (readProcessWithExitCode)

import Syntax (QCircuit(..), QGate(..), numQubitsQC)

--------------------------------------------------------------------------------
-- QCircuit -> OpenQASM 3

emit :: QCircuit -> IO String
emit qc = program (max 1 (numQubitsQC qc)) <$> emitGates qc

emitWithMeasurements :: QCircuit -> IO String
emitWithMeasurements qc =
  programWithMeasurements (max 1 (numQubitsQC qc)) <$> emitGates qc

emitGates :: QCircuit -> IO [String]
emitGates (QCircuit gates) = concat <$> mapM emitGate gates

emitGate :: QGate -> IO [String]
emitGate = \case
  QGPhase th            -> pure [gphase th]
  QH i                  -> pure [gh i]
  QS i                  -> pure [gs i]
  QX i                  -> pure [gx i]
  QCX c t               -> pure [gcx c t]
  QSWAP a b             -> pure [gswap a b]
  QP th i               -> pure [gp th i]
  QCP th c t            -> pure [gcp th c t]
  QDiagonalAt off n phis -> lowerDiagonal off n phis

--------------------------------------------------------------------------------
-- Gate string helpers

qref :: Natural -> String
qref i = "q[" ++ show i ++ "]"

gphase :: Double -> String
gphase th = "gphase(" ++ show th ++ ");"

gh, gs, gx :: Natural -> String
gh i = "h "    ++ qref i ++ ";"
gs i = "s "    ++ qref i ++ ";"
gx i = "x "    ++ qref i ++ ";"

gcx :: Natural -> Natural -> String
gcx c t = "cx " ++ qref c ++ ", " ++ qref t ++ ";"

gswap :: Natural -> Natural -> String
gswap a b = "swap " ++ qref a ++ ", " ++ qref b ++ ";"

gp :: Double -> Natural -> String
gp th i = "p(" ++ show th ++ ") " ++ qref i ++ ";"

gcp :: Double -> Natural -> Natural -> String
gcp th c t =
  "cp(" ++ show th ++ ") " ++ qref c ++ ", " ++ qref t ++ ";"

--------------------------------------------------------------------------------
-- Diagonal lowering via scripts.compiler.diag_to_qasm3
--
-- Protocol (JSON over stdin/stdout):
--   in:  {"off": int, "n": int, "phis": [float]}
--   out: {"gates": [{"name": str, "params": [float], "qubits": [int]}]}

data PyGate = PyGate
  { pgName   :: String
  , pgParams :: [Double]
  , pgQubits :: [Natural]
  } deriving (Eq, Show, Generic)

instance A.FromJSON PyGate where
  parseJSON = A.withObject "PyGate" $ \o ->
    PyGate <$> o A..: "name"
           <*> o A..: "params"
           <*> o A..: "qubits"

newtype PyProgram = PyProgram { pgGates :: [PyGate] }
  deriving (Eq, Show, Generic)

instance A.FromJSON PyProgram where
  parseJSON = A.withObject "PyProgram" $ \o ->
    PyProgram <$> o A..: "gates"

lowerDiagonal :: Natural -> Natural -> [Double] -> IO [String]
lowerDiagonal off n phis = do
  let payload = BL.unpack $ A.encode $ A.object
        [ "off"  A..= (fromIntegral off :: Integer)
        , "n"    A..= (fromIntegral n   :: Integer)
        , "phis" A..= phis
        ]
  (ec, out, err) <- readProcessWithExitCode
    "python3" ["-m", "scripts.compiler.diag_to_qasm3"] payload
  case ec of
    ExitFailure k ->
      error ("scripts.compiler.diag_to_qasm3 failed: " ++ show k ++ "\n" ++ err)
    ExitSuccess ->
      case A.eitherDecode (BL.pack out) of
        Left e  ->
          error ("scripts.compiler.diag_to_qasm3 bad JSON: " ++ e ++ "\n" ++ out)
        Right pr ->
          pure (concatMap renderGate (pgGates pr))

renderGate :: PyGate -> [String]
renderGate (PyGate nm ps qs) = case (nm, ps, qs) of
  ("gphase", [th], []   ) -> [gphase th]
  ("h",      [],   [i]  ) -> [gh    i]
  ("s",      [],   [i]  ) -> [gs    i]
  ("x",      [],   [i]  ) -> [gx    i]
  ("cx",     [],   [c,t]) -> [gcx   c t]
  ("swap",   [],   [a,b]) -> [gswap a b]
  ("p",      [th], [i]  ) -> [gp    th i]
  ("cp",     [th], [c,t]) -> [gcp   th c t]
  _ -> error
    ("unsupported gate from scripts.compiler.diag_to_qasm3: " ++ show (nm, ps, qs))

--------------------------------------------------------------------------------
-- Program templates

program :: Natural -> [String] -> String
program nq body = unlines $
  [ "OPENQASM 3;"
  , "include \"stdgates.inc\";"
  , "qubit[" ++ show nq ++ "] q;"
  , ""
  ] ++ body

programWithMeasurements :: Natural -> [String] -> String
programWithMeasurements nq body =
  program nq body ++ unlines
    ( [ "", "bit[" ++ show nq ++ "] c;" ]
      ++ [ "c[" ++ show i ++ "] = measure " ++ qref i ++ ";"
         | i <- [0 .. nq - 1] ]
    )

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
