import Prelude hiding ( zipWith, sum )
import qualified Prelude as P
import Data.Vector.Unboxed
import Data.Vector.Unboxed.Mutable

nNodes :: Int
nNodes = 5

nextNode :: Vector Int
nextNode = generate nNodes ((`mod` nNodes) . (+1))

xs, ys :: Vector Float
xs = generate nNodes ((*1.0) . fromIntegral)
ys = generate nNodes ((*0.0) . fromIntegral)

distance :: Int -> Int -> Float
distance i j = sqrt $ (xs!i - xs!j)^2 + (ys!i - ys!j)^2

totalDistance :: Vector Int -> Float
totalDistance v = sum $ zipWith distance v (backpermute v nextNode)

swapPair :: Int -> Int -> Vector Int -> Vector Int
swapPair i j v = modify (\v -> swap v i j) v

incCost :: Int -> Int -> Vector Int -> Float
incCost i j v = d + c - b - a
  where
    a = distance i (v!i)
    b = distance j (v!j)
    c = distance i j
    d = distance (v!i) (v!j)

expDv :: Int -> Int -> Vector Int -> Float -> Float -> Float -> Float
expDv i1 i2 v kB j t = exp(-j * (incCost i1 i2 v) / (kB * t))
