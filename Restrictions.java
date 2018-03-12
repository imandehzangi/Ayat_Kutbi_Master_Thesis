package rna_sequences;

import java.util.HashSet;
import java.util.Set;

public final class Restrictions {
	
	public int minMutationCount;
	
	public int maxMutationCount;
	
	public int minPairCount;
	
	public int maxPairCount;
	
	/** If <code>null</code>, all locations are valid. */
	public Set<Integer> validPairLocations;
	
	/** If <code>null</code>, all locations are valid. */
	public Set<Integer> validMutationLocations;
	
	public boolean useOrderedPairs;
	
	
	public Restrictions() {
		reset();
	}
	
	public void reset() {
		minMutationCount = 0;
		maxMutationCount = Integer.MAX_VALUE;
		minPairCount = 0;
		maxPairCount = Integer.MAX_VALUE;
		validMutationLocations = null;
		validPairLocations = null;
		useOrderedPairs = false;
	}
	
	public static Set<Integer> intSet(int... ns) {
		HashSet<Integer> set = new HashSet<>();
		for( int n : ns )
			set.add(n);
		return set;
	}
	
}
