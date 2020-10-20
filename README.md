# GABasedCombinatorialSearches
Code base to search optimal within NKp fitness labscape for optimal binary array
NOTE: The code requires python2.7, scientific python libraries numpy and scipy.  

1) Diversity.py: The code implements various flavors of genetic algorithm to find optimal solution of NKp landscape
2) generateCombination.py: The code implements functions to enumerate all possible
   binary arrays of a certain length that represent the combinatorial search space
3) MaxMin.py: The code implements maxmin selection algorithm to select diverse binary arrays to be used with in GA searches


Exaxmple:
To run the code: python Diversity.py

The variables of the NKp fitness landscape or GA searchs can be modified by the users within Diversity.py
