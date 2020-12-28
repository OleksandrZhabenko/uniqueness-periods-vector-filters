-- |
-- Module      :  Languages.UniquenessPeriods.Vector.Filters
-- Copyright   :  (c) OleksandrZhabenko 2020
-- License     :  MIT
-- Stability   :  Experimental
-- Maintainer  :  olexandr543@yahoo.com
--
-- A module allows to change the structure of the function output for the functions of
-- elements from 'RealFrac' class. At the moment only the equal intervals are supported.

module Languages.UniquenessPeriods.Vector.Filters (
  -- * One interval used
  intervalNRealFrac
  , unsafeTransfer1I5
  , transfer1IEq3
  -- * Several intervals
  , unsafeRearrangeIG
  , unsafeRearrangeIGV
  -- * Some basic usage examples
  , unsafeSwapIWithMaxI
  , unsafeSwapVecIWithMaxI
) where

import CaseBi
import qualified Data.Vector as V

-- | Given the minimum and maximum elements, a quantity of equal intervals, and an element in between the first two arguments (or equal to one of them), finds out the
-- index of the interval, to which the element belongs (starting from 1). The minimum element belongs to the interval with the index 1.
intervalNRealFrac
  :: RealFrac b => b
  -> b
  -> Int
  -> b
  -> Int
intervalNRealFrac minE maxE n x
 | maxE == minE = ceiling (0.5 * fromIntegral n)
 | otherwise = zero2One . ceiling $ fromIntegral n * (x - minE) / (maxE - minE)
{-# INLINE intervalNRealFrac #-}

zero2One :: Int -> Int
zero2One x = if x == 0 then 1 else x
{-# INLINE zero2One #-}

-- | Moves (if needed) the given value so that its result divides the new [min..max] interval in the same proportion as the starting one. Is intended to be used
-- for the arguments satisfying some additional constraints, but they are not checked (hence, its name prefix \"unsafe\"). For example, the second argument must be
-- greater than the first one, the fourth -- than the third one, and the fifth must be located in between the first two. Then the result is also located in between
-- the third and fourth arguments similarly.
unsafeTransfer1I5
  :: RealFrac b => b
  -> b
  -> b
  -> b
  -> b
  -> b
unsafeTransfer1I5 minE0 maxE0 minE1 maxE1 x
 | minE0 == maxE0 = x
 | otherwise = minE1 + (x - minE0) * (maxE1 - minE1) / (maxE0 - minE0)
{-# INLINE unsafeTransfer1I5 #-}

-- | A variant of the 'unsafeTransfer1I5' where the lengths of the both intervals (the old and the new ones) are equal.
transfer1IEq3
  :: RealFrac b => b
  -> b
  -> b
  -> b
transfer1IEq3 minE0 minE1 = (+ (minE1 - minE0))
{-# INLINE transfer1IEq3 #-}

-- | Makes a complex interval-based transformation moving the value from its own interval to the corresponding by the 'V.Vector' of tuples second element of the
-- respective pair with the first element being the starting number of the interval (numeration of them begins at 1). The 'V.Vector' argument must be sorted
-- by the first argument in the ascending order. Usually, its first elements in the tuples are from the range @[1..n]@. Number of the intervals are given as
-- the third argument and for many cases should not be greater than 10. There do exist several semantical constraints for the possible accurate arguments,
-- but they are not checked. For example, the first argument must be less than the second one; the fifth argument must be located between the first two ones;
-- the 'Int' argument must be greater than zero.
unsafeRearrangeIG
  :: RealFrac b => b
  -> b
  -> Int
  -> V.Vector (Int,Int)
  -> b
  -> b
unsafeRearrangeIG minE maxE n v x
 | minE == maxE = x
 | otherwise = x + fromIntegral (getBFst' (n0, v) n0 - n0) * (maxE - minE) / fromIntegral n
      where n0 = intervalNRealFrac minE maxE n x

-- | An unzipped variant of the 'unsafeRearrangeIG' function where the 'V.Vector' argument is internally 'V.zip'ped as the second argument with the 'V.Vector' @[1..n]@.
-- This allows to shorten the time of the arguments writing given only the resulting backpermutted indexes in the 'V.Vector'.
unsafeRearrangeIGV
  :: RealFrac b => b
  -> b
  -> Int
  -> V.Vector Int
  -> b
  -> b
unsafeRearrangeIGV minE maxE n v = unsafeRearrangeIG minE maxE n . V.zip (V.enumFromN 1 n) $ v
{-# INLINE unsafeRearrangeIGV #-}

-- | Swaps the k-th inner interval values with the maximum one's (that is the n-th one) values.
unsafeSwapIWithMaxI
  :: RealFrac b => b
  -> b
  -> Int -- ^ It is expected to be greater than 0, though this is not checked.
  -> Int -- ^ It is expected to be less than the previous argument, but greater than 0, though this is not checked.
  -> b -- ^ It is expected to lie between the first two arguments, though this is not checked.
  -> b
unsafeSwapIWithMaxI minE maxE n k = unsafeRearrangeIGV minE maxE n (V.generate n (\i -> if i == k - 1 then n - 1 else if i == n - 1 then k - 1 else i))
{-# INLINE unsafeSwapIWithMaxI #-}

-- | Swaps the inner intervals values (given by the 'V.Vector' of 'Int' elements that represent numbers-indices starting from 1 to n) with the maximum one's
-- (that is the n-th one) values. The 'V.Vector' must be not empty and sorted in the ascending order, though it is not checked. Be aware that this can
-- significantly change the density of the values and break some other properties for distributions.
unsafeSwapVecIWithMaxI
  :: RealFrac b => b
  -> b
  -> Int -- ^ It is expected to be greater than 0, though this is not checked.
  -> V.Vector Int -- ^ It is expected the 'V.Vector' to be sorted in the ascending order (indices are counted in it starting with 1 opposed to the usual behaviour for 'V.Vector's and are the numbers of the intervals in the range from 1 to n), and besides all the elements to be less than the previous argument, greater than 0 and to be not pairwise equal, though it is not checked.
  -> b -- ^ It is expected to lie between the first two arguments, though this is not checked.
  -> b
unsafeSwapVecIWithMaxI minE maxE n vI = unsafeRearrangeIGV minE maxE n (V.generate n h)
  where h i
         | getBFst' (False, V.zip (V.map (+ (-1)) vI) . V.replicate n $ True) i = n - 1
         | i == n - 1 = V.unsafeIndex vI 0
         | otherwise = i
{-# INLINE unsafeSwapVecIWithMaxI #-}
