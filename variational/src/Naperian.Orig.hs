{-# LINE 89 "apl.lhs" #-}
{- 

Code to accompany the paper
"APLicative Programming with Naperian Functors"
by Jeremy Gibbons
to appear in ESOP 2017

Requires GHC 7.10 or later

 -}

{-#  LANGUAGE StandaloneDeriving  #-}
{-#  LANGUAGE DeriveFunctor  #-}
{-#  LANGUAGE DeriveFoldable  #-}
{-#  LANGUAGE GADTs  #-}
{-#  LANGUAGE MultiParamTypeClasses  #-}
{-#  LANGUAGE FunctionalDependencies  #-}
{-#  LANGUAGE UndecidableInstances  #-}
{-#  LANGUAGE FlexibleInstances  #-}
{-#  LANGUAGE KindSignatures  #-}
{-#  LANGUAGE DataKinds  #-}
{-#  LANGUAGE PolyKinds  #-}
{-#  LANGUAGE TypeFamilies  #-}
{-#  LANGUAGE TypeOperators  #-}
{-#  LANGUAGE FlexibleContexts  #-}
{-#  LANGUAGE RankNTypes  #-}
{-#  LANGUAGE ScopedTypeVariables   #-}

import Prelude hiding (sum,lookup)
import Control.Applicative
import Data.Foldable (toList, fold)
import Data.Monoid (Sum(..),Product(..))
import Data.List (intersperse)
import GHC.TypeLits (Symbol)
import Data.Array (Array, listArray, accumArray, ixmap, bounds, elems, Ix(..))
{-# LINE 161 "apl.lhs" #-}

{-# LINE 578 "apl.lhs" #-}
data Nat :: *   where
  Z  ::         Nat
  S  :: Nat ->  Nat
{-# LINE 603 "apl.lhs" #-}
type One    = 'S 'Z
type Two    = 'S One
type Three  = 'S Two
type Four   = 'S Three
{-# LINE 611 "apl.lhs" #-}
data Vector :: Nat -> * -> *   where
  VNil   ::                     Vector 'Z a
  VCons  :: a -> Vector n a ->  Vector ('S n) a
{-# LINE 620 "apl.lhs" #-}
v123 :: Vector Three Int
v123  = VCons 1 (VCons 2 (VCons 3 VNil))
{-# LINE 638 "apl.lhs" #-}
vmap :: (a->b) -> Vector n a -> Vector n b
vmap f VNil          = VNil
vmap f (VCons x xs)  = VCons (f x) (vmap f xs)
{-# LINE 644 "apl.lhs" #-}
v456 :: Vector Three Int
v456  = vmap (\ x -> 3+x) v123
{-# LINE 651 "apl.lhs" #-}
vtail :: Vector (S n) a -> Vector n a
vtail (VCons x xs) = xs
{-# LINE 659 "apl.lhs" #-}
vzipWith :: (a->b->c) -> Vector n a -> Vector n b -> Vector n c
vzipWith f VNil         VNil         = VNil
vzipWith f (VCons a x)  (VCons b y)  = VCons (f a b) (vzipWith f x y)
{-# LINE 695 "apl.lhs" #-}
class Count (n::Nat) where 
  vreplicate  :: a -> Vector n a
  vlength     :: Vector n a -> Int
{-# LINE 703 "apl.lhs" #-}
instance Count 'Z where 
  vreplicate a          = VNil
  vlength xs            = 0

instance Count n => Count ('S n) where 
  vreplicate a          = VCons a (vreplicate a)
  vlength xs            = 1 + vlength (vtail xs)
{-# LINE 741 "apl.lhs" #-}
instance (Show a) => Show (Vector n a) where
  show = glue . map show . toList

glue :: [String] -> String
glue x = "<" ++ concat (intersperse "," x) ++ ">"
{-# LINE 766 "apl.lhs" #-}
data Pair a = P a a
{-# LINE 799 "apl.lhs" #-}
areplicate :: Applicative f => a -> f a
areplicate = pure
{-# LINE 812 "apl.lhs" #-}
azipWith :: Applicative f => (a -> b -> c) -> f a -> f b -> f c
azipWith h xs ys = (pure h <*> xs) <*> ys
{-# LINE 818 "apl.lhs" #-}
instance Functor (Vector n) where 
  fmap = vmap

instance Count n => Applicative (Vector n) where
  pure   = vreplicate 
  (<*>)  = vzipWith (\ f x -> f x)
{-# LINE 833 "apl.lhs" #-}
instance Functor Pair where
  fmap f (P x y) = P (f x) (f y)

instance Applicative Pair where
  pure x           = P x x
  P f g <*> P x y  = P (f x) (g y)
{-# LINE 893 "apl.lhs" #-}
class Functor f => Naperian f where
  type Log f
  lookup      :: f a -> (Log f -> a)
  tabulate    :: (Log f -> a) -> f a
  positions   :: f (Log f)

  tabulate h  = fmap h positions
  positions   = tabulate id
{-# LINE 955 "apl.lhs" #-}
transpose :: (Naperian f, Naperian g) => f (g a) -> g (f a)
transpose = tabulate . fmap tabulate . flip . fmap lookup . lookup
{-# LINE 980 "apl.lhs" #-}
transposeF :: (Functor f, Naperian g) => f (g a) -> g (f a)
transposeF x = tabulate (\ p -> fmap ($ p) (fmap lookup x))
{-# LINE 992 "apl.lhs" #-}
instance Naperian Pair where 
  type Log Pair     = Bool
  lookup (P x y) b  = if b then y else x
  positions         = P False True
{-# LINE 1001 "apl.lhs" #-}
data Fin :: Nat -> *   where
  FZ ::           Fin ('S n)
  FS :: Fin n ->  Fin ('S n)
{-# LINE 1013 "apl.lhs" #-}
vlookup :: Vector n a -> Fin n -> a
vlookup (VCons a x) FZ      = a
vlookup (VCons a x) (FS n)  = vlookup x n
{-# LINE 1022 "apl.lhs" #-}
viota :: Count n => Vector n (Fin n) 
viota = viota' (vreplicate ()) where
  viota' :: Vector m () -> Vector m (Fin m)
  viota' VNil           = VNil
  viota' (VCons () xs)  = VCons FZ (fmap FS (viota' xs))
{-# LINE 1030 "apl.lhs" #-}
instance Count n => Naperian (Vector n) where
  type Log (Vector n)  = Fin n
  lookup               = vlookup
  positions            = viota 
{-# LINE 1041 "apl.lhs" #-}
twoXthree :: Vector Two (Vector Three Int)
twoXthree = fmap (fmap succ) (VCons v123 (VCons v456 VNil))
{-# LINE 1064 "apl.lhs" #-}
instance Foldable (Vector n) where
  foldr f e VNil          = e
  foldr f e (VCons x xs)  = f x (foldr f e xs)
{-# LINE 1070 "apl.lhs" #-}
instance Foldable Pair where
  foldr f e (P x y) = f x (f y e)
{-# LINE 1082 "apl.lhs" #-}
sum :: (Num a, Foldable t) => t a -> a
sum = foldr (+) 0
{-# LINE 1104 "apl.lhs" #-}
data State s a = State { runState :: s -> (a,s) }
{-# LINE 1113 "apl.lhs" #-}
instance Functor (State s) where
  fmap f (State h) = State ( \ s -> let (x,s') = h s in (f x, s'))

instance Applicative (State s) where
  pure a               = State (\ s -> (a,s))
  State h <*> State k  = State (\ s -> let (f,s') = h s ; (x,s'') = k s' in (f x, s''))
{-# LINE 1127 "apl.lhs" #-}
increase :: Num a => a -> State a a
increase m = State (\ n -> (m+n,m+n))
{-# LINE 1134 "apl.lhs" #-}
sums :: (Num a, Traversable t) => t a -> t a
sums xs = fst (runState (traverse increase xs) 0)
{-# LINE 1148 "apl.lhs" #-}
instance Traversable Pair where
  traverse f (P x y) = (pure P <*> f x) <*> f y

instance Traversable (Vector n) where
  traverse f VNil          = pure VNil
  traverse f (VCons x xs)  = (pure VCons <*> f x) <*> traverse f xs
{-# LINE 1157 "apl.lhs" #-}
class (Applicative f, Naperian f, Traversable f) => Dimension f where
  size :: f a -> Int
  size = length . toList
{-# LINE 1165 "apl.lhs" #-}
instance             Dimension Pair        where size = const 2
instance Count n =>  Dimension (Vector n)  where size = vlength
{-# LINE 1171 "apl.lhs" #-}
data Perfect :: Nat -> * -> *   where
  Leaf  :: a ->                   Perfect 'Z a
  Bin   :: Pair (Perfect n a) ->  Perfect ('S n) a  
{-# LINE 1179 "apl.lhs" #-}
data Binary :: *   where
  Unit     ::            Binary
  Twice    :: Binary ->  Binary
  TwiceP1  :: Binary ->  Binary
{-# LINE 1186 "apl.lhs" #-}
bin2int :: Binary -> Int
bin2int Unit         = 1
bin2int (Twice n)    = 2 * bin2int n
bin2int (TwiceP1 n)  = 2 * bin2int n + 1
{-# LINE 1194 "apl.lhs" #-}
data BVector :: Binary -> * -> *   where
  VSingle  :: a ->                                 BVector 'Unit a
  VJoin    ::       BVector n a -> BVector n a ->  BVector ('Twice n) a
  VJoinP1  :: a ->  BVector n a -> BVector n a ->  BVector ('TwiceP1 n) a
{-# LINE 1223 "apl.lhs" #-}
innerp :: (Num a, Dimension f) => f a -> f a -> a
innerp xs ys = sum (azipWith (*) xs ys)
{-# LINE 1230 "apl.lhs" #-}
matrixp ::  (Num a, Dimension f, Dimension g, Dimension h) =>
            f (g a) -> g (h a) -> f (h a)
matrixp xss yss = azipWith (azipWith innerp)  (fmap areplicate xss) 
                                              (areplicate (transpose yss))
{-# LINE 1258 "apl.lhs" #-}
vv123456 :: Vector Two (Vector Three Int)
vv123456 =  VCons v123 (VCons v456 VNil)
{-# LINE 1276 "apl.lhs" #-}
data Hyper0 :: * -> *   where  --  to be refined later
  Scalar0 ::              a ->                    Hyper0 a
  Prism0  ::  Count n =>  Hyper0 (Vector n a) ->  Hyper0 a
{-# LINE 1309 "apl.lhs" #-}
data Hyper1 :: Nat -> * -> *   where  --  to be refined later
  Scalar1 ::              a ->                      Hyper1 'Z a
  Prism1  ::  Count n =>  Hyper1 r (Vector n a) ->  Hyper1 ('S r) a
{-# LINE 1331 "apl.lhs" #-}
data Hyper2 :: [Nat] -> * -> *   where  --  to be refined later
  Scalar2 ::              a ->                       Hyper2 '[] a
  Prism2  ::  Count n =>  Hyper2 ns (Vector n a) ->  Hyper2 (n ': ns) a
{-# LINE 1362 "apl.lhs" #-}
class Shapely fs where
  hreplicate  :: a -> Hyper fs a
  hsize       :: Hyper fs a -> Int
{-# LINE 1368 "apl.lhs" #-}
instance Shapely '[]  where
  hreplicate a      = Scalar a
  hsize             = const 1
instance forall f fs . (Dimension f, Shapely fs) => Shapely (f ': fs) where
  hreplicate a      = Prism (hreplicate (areplicate a))
  hsize (Prism x)   = size (first x) * hsize x
{-# LINE 1380 "apl.lhs" #-}
first :: Shapely fs => Hyper fs a -> a
first (Scalar a)  = a
first (Prism x)   = head (toList (first x))
{-# LINE 1392 "apl.lhs" #-}
data Hyper :: [* -> *] -> * -> *   where  --  final version
  Scalar ::                                a ->               Hyper '[] a
  Prism  ::  (Dimension f, Shapely fs) =>  Hyper fs (f a) ->  Hyper (f ': fs) a
{-# LINE 1398 "apl.lhs" #-}
h123456 :: Hyper '[Vector Three, Vector Two] Int
h123456 = Prism (Prism (Scalar vv123456))
{-# LINE 1403 "apl.lhs" #-}
instance Functor (Hyper fs) where
  fmap f (Scalar a)  = Scalar (f a)
  fmap f (Prism x)   = Prism (fmap (fmap f) x)
{-# LINE 1411 "apl.lhs" #-}
hzipWith :: (a->b->c) -> Hyper fs a -> Hyper fs b -> Hyper fs c
hzipWith f (Scalar a)  (Scalar b)  = Scalar (f a b)
hzipWith f (Prism x)   (Prism y)   = Prism (hzipWith (azipWith f) x y)
{-# LINE 1417 "apl.lhs" #-}
instance Shapely fs => Applicative (Hyper fs) where
  pure   = hreplicate 
  (<*>)  = hzipWith (\ f x -> f x)
{-# LINE 1428 "apl.lhs" #-}
reduceBy :: (a->a->a, a) -> Hyper (f ': fs) a -> Hyper fs a
reduceBy (f,e) (Prism x) = fmap (foldr f e) x
{-# LINE 1434 "apl.lhs" #-}
transposeHyper :: Hyper (f ': (g ': fs)) a -> Hyper (g ': (f ': fs)) a
transposeHyper (Prism (Prism x)) = Prism (Prism (fmap transpose x)) 
{-# LINE 1448 "apl.lhs" #-}
instance Show a => Show (Hyper fs a) where
  show = showHyper . fmap show
    where
      showHyper :: Hyper gs String -> String
      showHyper (Scalar s) = s
      showHyper (Prism x) = showHyper (fmap (glue . toList) x)
{-# LINE 1457 "apl.lhs" #-}
t6a :: Hyper '[] Int
t6a = Scalar 3
t6b :: Hyper ((Vector Three) ': '[]) Int
t6b = Prism (Scalar v123)
t6c :: Hyper '[] (Vector Two (Vector Three Int))
t6c = Scalar twoXthree
t6d :: Hyper ((Vector Two) ': '[]) (Vector Three Int)
t6d = Prism (Scalar twoXthree)
t6e :: Hyper ((Vector Three) ': ((Vector Two) ': '[])) Int
t6e = Prism (Prism (Scalar twoXthree))
t6f :: Hyper ((Vector Four) ': ((Vector Three) ': ((Vector Two) ': '[]))) Int
t6f = Prism (fmap v t6e) where
  v :: Int -> Vector Four Int
  v n = VCons m (VCons (m+1) (VCons (m+2) (VCons (m+3) VNil))) where m = 4*n-7
t6g :: Hyper ((Vector Four) ': '[]) Int
t6g = Prism (Scalar (VCons 1 $ VCons 2 $ VCons 3 $ VCons 4 VNil))
{-# LINE 1482 "apl.lhs" #-}
unary :: Shapely fs => (a->b) -> (Hyper fs a -> Hyper fs b)
unary = fmap
{-# LINE 1492 "apl.lhs" #-}
class (Shapely fs, Shapely gs) => Alignable fs gs where
  align :: Hyper fs a -> Hyper gs a
{-# LINE 1497 "apl.lhs" #-}
instance Alignable '[] '[] where
  align = id
{-# LINE 1502 "apl.lhs" #-}
instance (Dimension f, Alignable fs gs) => Alignable (f ': fs) (f ': gs) where
  align (Prism x) = Prism (align x)
{-# LINE 1509 "apl.lhs" #-}
instance (Dimension f, Shapely fs) => Alignable '[] (f ': fs) where
  align (Scalar a) = hreplicate a
{-# LINE 1539 "apl.lhs" #-}
type family Max (fs :: [* -> *]) (gs :: [* -> *]) :: [* -> *] where
  Max '[]            '[]            = '[]
  Max '[]            (f ': gs)  = (f ': gs)
  Max (f ': fs)  '[]            = (f ': fs)
  Max (f ': fs)  (f ': gs)  = (f ': Max fs gs)
{-# LINE 1572 "apl.lhs" #-}
binary0 :: --  to be refined later
  (Max fs gs ~ hs, Alignable fs hs, Alignable gs hs) => 
  (a->b->c) -> (Hyper fs a -> Hyper gs b -> Hyper hs c)
binary0 f x y = hzipWith f (align x) (align y)
{-# LINE 1603 "apl.lhs" #-}
type family IsCompatible (fs :: [* -> *]) (gs :: [* -> *]) :: IsDefined Symbol where
  IsCompatible '[]            '[]            = Defined
  IsCompatible '[]            (f ': gs)  = Defined
  IsCompatible (f ': fs)  '[]            = Defined
  IsCompatible (f ': fs)  (f ': gs)  = IsCompatible fs gs
  IsCompatible (f ': fs)  (g ': gs)  = Undefined "Mismatching dimensions"
{-# LINE 1615 "apl.lhs" #-}
data IsDefined e = Defined | Undefined e
{-# LINE 1619 "apl.lhs" #-}
binary :: --  final version
  (  IsCompatible fs gs ~ Defined, Max fs gs ~ hs, Alignable fs hs, Alignable gs hs) => 
  (a->b->c) -> (Hyper fs a -> Hyper gs b -> Hyper hs c)
binary f x y = binary0 f x y
{-# LINE 1637 "apl.lhs" #-}
t7b :: Hyper '[Vector Two] Int
t7b = reduceBy ((+),0) t6e
t7c :: Hyper '[Vector Three] Int
t7c = reduceBy ((+),0) (transposeHyper t6e)
{-# LINE 1644 "apl.lhs" #-}
matrixp0 ::  (Num a, Dimension f, Dimension g, Dimension h) =>
            Hyper '[ g,f] a -> Hyper '[ h,g] a -> Hyper '[ h,f] a
matrixp0 xss yss = case (xss, transposeHyper yss) of (Prism xss', Prism (Prism yss')) -> hzipWith innerp (Prism (fmap areplicate xss')) (Prism (Prism (fmap areplicate yss')))
{-# LINE 1665 "apl.lhs" #-}
data HyperR :: [* -> *] -> * -> *   where
  ScalarR  ::                                a ->                HyperR '[] a
  PrismR   ::  (Dimension f, Shapely fs) =>  HyperR fs (f a) ->  HyperR (f ': fs) a
  ReplR    ::  (Dimension f, Shapely fs) =>  HyperR fs a ->      HyperR (f ': fs) a
{-# LINE 1681 "apl.lhs" #-}
instance Functor (HyperR fs) where
  fmap f (ScalarR a)  = ScalarR (f a)
  fmap f (PrismR x)   = PrismR (fmap (fmap f) x)
  fmap f (ReplR x)    = ReplR (fmap f x)
{-# LINE 1691 "apl.lhs" #-}
rzipWith :: Shapely fs => (a->b->c) -> HyperR fs a -> HyperR fs b -> HyperR fs c
rzipWith f (ScalarR a) (ScalarR b)  = ScalarR (f a b)
rzipWith f (PrismR x) (PrismR y)    = PrismR (rzipWith (azipWith f) x y)
rzipWith f (PrismR x) (ReplR y)     = PrismR (rzipWith (azipWithL f) x y)
rzipWith f (ReplR x) (PrismR y)     = PrismR (rzipWith (azipWithR f) x y)
rzipWith f (ReplR x) (ReplR y)      = ReplR (rzipWith f x y)
{-# LINE 1701 "apl.lhs" #-}
azipWithL :: Functor f => (a->b->c) -> f a -> b -> f c
azipWithL f xs y = fmap (\ x -> f x y) xs

azipWithR :: Functor f => (a->b->c) -> a -> f b -> f c
azipWithR f x ys = fmap (\ y -> f x y) ys
{-# LINE 1715 "apl.lhs" #-}
rtranspose ::  (Shapely fs, Dimension f, Dimension g) => 
               HyperR (f ': g ': fs) a ->  HyperR (g ': f ': fs) a
rtranspose (PrismR (PrismR x))  = PrismR (PrismR (fmap transpose x))
rtranspose (PrismR (ReplR x))   = ReplR (PrismR x)
rtranspose (ReplR (PrismR x))   = PrismR (ReplR x)
rtranspose (ReplR (ReplR x))    = ReplR (ReplR x)
{-# LINE 1732 "apl.lhs" #-}
forceReplR :: Shapely fs => HyperR fs a -> HyperR fs a
forceReplR (ReplR x)  = PrismR (fmap areplicate x)
forceReplR x          = x
{-# LINE 1739 "apl.lhs" #-}
data HyperT :: [* -> *] -> * -> *   where
  ScalarT  ::  a ->                                        HyperT '[] a
  PrismT   ::  (Dimension f, Shapely fs) =>  
               HyperT fs (f a) ->                          HyperT (f ': fs) a
  TransT   ::  (Dimension f, Dimension g, Shapely fs) => 
               HyperT (f ': g ': fs) a ->  HyperT (g ': f ': fs) a
{-# LINE 1752 "apl.lhs" #-}
transT   ::  (Dimension f, Dimension g, Shapely fs) => 
             HyperT (f ': g ': fs) a ->           HyperT (g ': f ': fs) a
transT (TransT x) = x
transT x = TransT x
{-# LINE 1768 "apl.lhs" #-}
forceTransT ::  (Dimension f, Dimension g, Shapely fs) => 
                HyperT (f ': g ': fs) a ->           HyperT (f ': g ': fs) a
forceTransT (TransT (PrismT (PrismT x)))
               =  PrismT (PrismT (fmap transpose x))
forceTransT (TransT (PrismT x@(TransT _)))
               =  case forceTransT x of 
                    PrismT x' -> PrismT (PrismT (fmap transpose x'))
forceTransT x  =  x
{-# LINE 1785 "apl.lhs" #-}
tzipWith ::  Shapely fs => 
             (a->b->c) -> HyperT fs a -> HyperT fs b -> HyperT fs c
tzipWith f (ScalarT a) (ScalarT b)  = ScalarT (f a b)
tzipWith f (PrismT x) (PrismT y)    = PrismT (tzipWith (azipWith f) x y)
tzipWith f (TransT x) (TransT y)    = TransT (tzipWith f x y) 
tzipWith f x@(TransT _) (PrismT y)  = tzipWith f (forceTransT x) (PrismT y)
tzipWith f (PrismT x) y@(TransT _)  = tzipWith f (PrismT x) (forceTransT y)
{-# LINE 1798 "apl.lhs" #-}
deriving instance Functor (HyperT fs)
{-# LINE 1812 "apl.lhs" #-}
data HyperS :: [* -> *] -> * -> * where
  ScalarS  ::  a ->                                        HyperS '[] a
  PrismS   ::  (Dimension f, Shapely fs) =>  
               HyperS fs (f a) ->                          HyperS (f ': fs) a
  ReplS    ::  (Dimension f, Shapely fs) =>  
               HyperS fs a ->                              HyperS (f ': fs) a
  TransS   ::  (Dimension f, Dimension g, Shapely fs) => 
               HyperS (f ': g ': fs) a ->  HyperS (g ': f ': fs) a
{-# LINE 1823 "apl.lhs" #-}
szipWith ::  Shapely fs => 
             (a->b->c) -> HyperS fs a -> HyperS fs b -> HyperS fs c
szipWith f  (ScalarS a) (ScalarS b)  = ScalarS (f a b)
szipWith f  (PrismS x) (PrismS y)    = PrismS (szipWith (azipWith f) x y)
szipWith f  (PrismS x)  (ReplS y)    = PrismS (szipWith (azipWithL f) x y)
szipWith f  (ReplS x)  (PrismS y)    = PrismS (szipWith (azipWithR f) x y)
szipWith f  (ReplS x)  (ReplS y)     = ReplS (szipWith f x y)
szipWith f  (TransS x) (TransS y)    = TransS (szipWith f x y)
szipWith f  (TransS x) 
            (PrismS (ReplS y))       = TransS (szipWith f x (ReplS (PrismS y)))
szipWith f  (TransS x) 
            (ReplS (PrismS y))       = TransS (szipWith f x (PrismS (ReplS y)))
szipWith f  (TransS x) 
            (ReplS (ReplS y))        = TransS (szipWith f x (ReplS (ReplS y)))
szipWith f  x@(TransS _) 
            y@(PrismS (PrismS _))    = szipWith f (forceTransS x) y
szipWith f  x@(TransS _) 
            y@(PrismS (TransS _))    = szipWith f (forceTransS x) y
szipWith f  x@(TransS _) 
            y@(ReplS (TransS _))     = szipWith f (forceTransS x) y
szipWith f  x y@(TransS _)           = szipWith (flip f) y x
{-# LINE 1847 "apl.lhs" #-}
forceTransS :: (Dimension f, Dimension g, Shapely fs) => 
               HyperS (f ': g ': fs) a ->           HyperS (f ': g ': fs) a
{-# LINE 1856 "apl.lhs" #-}
forceReplS :: Shapely fs => HyperS fs a -> HyperS fs a
forceReplS (ReplS x) = PrismS (fmap areplicate x)
forceReplS x = x

forceTransS (TransS (ReplS (ReplS x))) = ReplS (ReplS x)
forceTransS (TransS (PrismS (ReplS x))) = ReplS (PrismS x)
forceTransS (TransS (ReplS (PrismS x))) =  PrismS (ReplS x)
forceTransS (TransS (PrismS (PrismS x))) = PrismS (PrismS (fmap transpose x))
forceTransS (TransS (PrismS x@(TransS _))) = forceTransS (TransS (PrismS (forceTransS x)))
forceTransS (TransS (ReplS x@(TransS _))) = forceTransS (TransS (ReplS (forceTransS x)))
forceTransS x = x

instance Shapely fs => Functor (HyperS fs) where
  fmap f (ScalarS a)  = ScalarS (f a)
  fmap f (PrismS x)   = PrismS (fmap (fmap f) x)
  fmap f (ReplS x)    = ReplS (fmap f x)
  fmap f (TransS x)   = TransS (fmap f x)

instance Shapely fs => Foldable (HyperS fs) where
  foldr f e (ScalarS a)    = f a e
  foldr f e (PrismS x)     = foldr (flip (foldr f)) e x
  foldr f e x@(ReplS _)    = foldr f e (forceReplS x)
  foldr f e x@(TransS _)   = foldr f e (forceTransS x)

instance Shapely fs => Traversable (HyperS fs) where
  traverse f (ScalarS a)          = pure ScalarS <*> f a
  traverse f (PrismS x)  = pure PrismS <*> traverse (traverse f) x
  traverse f x@(ReplS _)  = traverse f (forceReplS x)
  traverse f x@(TransS _)  = traverse f (forceTransS x)

transS   ::  (Dimension f, Dimension g, Shapely fs) => 
               HyperS (f ': g ': fs) a ->           HyperS (g ': f ': fs) a
transS (TransS x) = x
transS x = TransS x

hStranspose   ::  (Dimension f, Dimension g, Shapely fs) => 
               HyperS (f ': g ': fs) a ->           HyperS (g ': f ': fs) a
hStranspose = transS
{-# LINE 1927 "apl.lhs" #-}
elements :: Shapely fs => Hyper fs a -> [a]
elements (Scalar a)  = [a]
elements (Prism x)   = concat (map toList (elements x))
{-# LINE 1939 "apl.lhs" #-}
data Flat fs a where 
  Flat :: Shapely fs => Array Int a -> Flat fs a
{-# LINE 1945 "apl.lhs" #-}
flatten :: Shapely fs => Hyper fs a -> Flat fs a
flatten x = Flat (listArray (0, hsize x - 1) (elements x))
{-# LINE 1971 "apl.lhs" #-}
data Sparse fs a where 
  Sparse :: Shapely fs => a -> Array Int (Int,a) -> Sparse fs a
{-# LINE 1979 "apl.lhs" #-}
unsparse :: forall fs a . Shapely fs => Sparse fs a -> Flat fs a
unsparse x@(Sparse e xs) = Flat (accumArray second e (0, l -1) (elems xs))
  where  l = hsize (hreplicate () :: Hyper fs ())
         second b a = a
{-# LINE 2023 "apl.lhs" #-}
t8a :: Flat '[Vector Three, Vector Two] Int
t8a = flatten t6e
{-# LINE 2061 "apl.lhs" #-}
type family Add (m :: Nat) (n :: Nat) :: Nat where
  Add 'Z n      = n
  Add ('S m) n  = 'S (Add m n)

vappend :: Vector m a -> Vector n a -> Vector (Add m n) a
vappend VNil y         = y
vappend (VCons a x) y  = VCons a (vappend x y)

vsplit :: Count m => Vector (Add m n) a -> (Vector m a, Vector n a)
vsplit = vsplit' (vreplicate ()) where
  vsplit' :: Vector m () -> Vector (Add m n) a -> (Vector m a, Vector n a)
  vsplit' VNil          x            = (VNil, x)
  vsplit' (VCons () m)  (VCons a x)  = let (y,z) = vsplit' (m) x in (VCons a y, z)

type family Mult (m :: Nat) (n :: Nat) :: Nat where
  Mult 'Z n      = 'Z
  Mult ('S m) n  = Add n (Mult m n)

vconcat :: Vector m (Vector n a) -> Vector (Mult m n) a
vconcat VNil          = VNil
vconcat (VCons x xs)  = vappend x (vconcat xs)

vgroup :: (Count m, Count n) => Vector (Mult m n) a -> Vector m (Vector n a)
vgroup = vgroup' (vreplicate ()) where
  vgroup' :: Count n => Vector m () -> Vector (Mult m n) a -> Vector m (Vector n a)
  vgroup' VNil          VNil  = VNil
  vgroup' (VCons () m)  x     = let (y,z) = vsplit x in VCons y (vgroup' m z)
