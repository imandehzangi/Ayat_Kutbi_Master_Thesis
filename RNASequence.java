package rna_sequences;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Predicate;

// immutable
public class RNASequence {
	
	private Node[] nodes;
	private Set<Pair> pairs;
	
	
	RNASequence(Letter[] letters, Collection<Integer> mutations, Set<Pair> pairs) {
		nodes = new Node[letters.length];
		this.pairs = Collections.unmodifiableSet(pairs);
		
		for( Integer i : mutations )
			nodes[i] = new Node(this, i, letters[i], true);
		for( int i = 0; i < nodes.length; i++ )
			if( nodes[i] == null )
				nodes[i] = new Node(this, i, letters[i], false);
		
		for( Pair p : pairs ) {
			nodes[p.getStart()].setPair(p);
			nodes[p.getEnd()].setPair(p);
		}
	}
	
	private static <T> Collection<Set<T>> subsets(List<T> list, int minCardinality, int maxCardinality, Predicate<Set<T>> filter) {
		if( maxCardinality > list.size() )
			maxCardinality = list.size();
		
		ArrayList<Set<T>> subsets = new ArrayList<>();
		for( int size = minCardinality; size <= maxCardinality; size++ ) {
			// create all subsets of size 'size'
			// -- Initialize iterators
			int[] i = new int[size];
			for( int a = 0; a < i.length; a++ )
				i[a] = a;
			
			while(true) {
				// -- build set
				Set<T> subset = new HashSet<>();
				for( int index : i )
					subset.add( list.get(index) );
				
				// -- validate set
				if( filter.test(subset) )
					subsets.add(subset);
				
				// -- iterate 'i'
				int a;
				for( a = size - 1; a >= 0; a-- ) 
					if( ++i[a] < list.size() - (i.length - 1 - a) )
						break;
				if( a < 0 )
					break; // break while
				for( int b = a + 1; b < size; b++ )
					i[b] = i[b - 1] + 1;
			}
		}
		return subsets;
	}
	
	private static <T> Collection<Set<T>> subsets(List<T> list, int minCardinality, int maxCardinality) {
		return subsets(list, minCardinality, maxCardinality, (set) -> {
			return true;
		});
	}
	
	/**
	 * Generate all of the possible sequence given the input nucleotide sequence
	 * and restrictions.
	 * 
	 * @param letters
	 *            A nucleotide sequence
	 * @param restrictions
	 *            Used to restrict the possible output sequence
	 * @return A collection of all of the RNA Sequences that satisfy the given
	 *         restrictions over the given nucleotide sequence. The elements in
	 *         the collection will not be in any particular order, but iterating
	 *         over the sequence will always happen in the same order.
	 */
	public static Collection<RNASequence> possibleSequences(Letter[] letters, Restrictions restrictions) {
		int length = letters.length;
		
		int minMutationCount = Math.max( 0, restrictions.minMutationCount );
		int maxMutationCount = Math.min( length, restrictions.maxMutationCount );
		int minPairCount = Math.max( 0, restrictions.minPairCount );
		int maxPairCount = Math.min( length / 2, restrictions.maxPairCount );
		Set<Integer> validPairLocations = restrictions.validPairLocations;
		Set<Integer> validMutationLocations = restrictions.validMutationLocations;
		
		if( validPairLocations == null ) {
			validPairLocations = new HashSet<>();
			for( int i = 0; i < length; i++ )
				validPairLocations.add(i);
		}
		if( validMutationLocations == null ) {
			validMutationLocations = new HashSet<>();
			for( int i = 0; i < length; i++ )
				validMutationLocations.add(i);
		}
		
		// calculate all the possible mutation sets of the given length
		ArrayList<Integer> possibleMutationLocations = new ArrayList<>(length);
		for( int i = 0; i < length; i++ )
			if( validMutationLocations.contains(i) )
				possibleMutationLocations.add(i);
		Collection<Set<Integer>> mutationSets = subsets(possibleMutationLocations, minMutationCount, maxMutationCount);
		
		LinkedList<RNASequence> seqs = new LinkedList<>();
		for( Set<Integer> mutationSet : mutationSets ) {
			if( mutationSet.size() == 0 ) {
				// BASE CASE (for recursion): don't split the sequence
				
				// calculate the set of all possible pairing locations
				ArrayList<Pair> possiblePairs = new ArrayList<>(length * length - length);
				for( int n1 = 0; n1 < length; n1++ )
					for( int n2 = n1 + 1; n2 < length; n2++ )
						if( validPairLocations.contains(n1) && validPairLocations.contains(n2) && letters[n1].isCompatible(letters[n2]) ) {
							possiblePairs.add( new Pair(n1, n2) );
							if( restrictions.useOrderedPairs )
								possiblePairs.add( new Pair(n2, n1) );
						}
				
				// calculate all the valid subsets of possible pairs (ignoring mutations)
				Collection<Set<Pair>> pairSets = subsets(possiblePairs, minPairCount, maxPairCount, (set) -> {
					boolean[] flags = new boolean[length];
					for( Pair pair : set ) {
						if( flags[pair.getStart()] || flags[pair.getEnd()] )
							return false;
						flags[pair.getStart()] = true;
						flags[pair.getEnd()] = true;
					}
					return true;
				});
				
				// add each pair,mutation set combination to the collection of valid sequences
				for( Set<Pair> pairSet : pairSets )
					seqs.add( new RNASequence(letters, mutationSet, pairSet) );
				
			} else {  // mutationSet.size() != 0
				// split the sequence at each mutation since binds can not cross over mutations
				
				ArrayList<Integer> ends = new ArrayList<>( mutationSet.size() );
				ends.addAll(mutationSet);
				ends.add(length);
				
				TreeMap<Integer, ArrayList<RNASequence>> subSeqs = new TreeMap<>();
				
				{	int start = 0;
					for( int end : ends ) {
						// grab a sub-sequence between start (inclusive) and end (exclusive), and recurse
						Letter[] subLetter = new Letter[end - start];
						System.arraycopy(letters, start, subLetter, 0, subLetter.length);
	
						Restrictions subRestrictions = new Restrictions();
						subRestrictions.minMutationCount = 0;
						subRestrictions.maxMutationCount = 0;
						subRestrictions.minPairCount = 0;
						subRestrictions.maxPairCount = maxPairCount;
						subRestrictions.validMutationLocations = new HashSet<>();
						subRestrictions.validPairLocations  = new HashSet<>();
						for( int loc : validPairLocations ) {
							int newLoc = loc - start;
							if( 0 <= newLoc && newLoc < subLetter.length )
								subRestrictions.validPairLocations.add(newLoc);
						}
						subRestrictions.useOrderedPairs = restrictions.useOrderedPairs;
	
						Collection<RNASequence> c = possibleSequences(subLetter, subRestrictions);
						subSeqs.put( start, new ArrayList<>(c) );
	
						// iterate
						start = end + 1;
					}
				}
				
				// for all valid combinations of 'subSeqs', join them into one RNASequence and add to the collection of valid sequence, 'seqs'
				ArrayList<Integer> starts = new ArrayList<>( subSeqs.keySet() );
				int[] i = new int[starts.size()];
				while(true) {
					// create RNASequence if valid, and add it to 'seqs'; otherwise, skip
					Set<Pair> pairSet = new HashSet<>();
					for( int k = 0; k < i.length; k++ ) {
						int start = starts.get(k);
						Set<Pair> subPairs = subSeqs.get(start).get(i[k]).getPairs();
						for( Pair p : subPairs )
							pairSet.add( p.shift(start) );
					}
					
					if( minPairCount <= pairSet.size() && pairSet.size() <= maxPairCount )
						seqs.add( new RNASequence(letters, mutationSet, pairSet) );
					
					// iterate
					int k;
					for( k = i.length - 1; k >= 0; k-- ) {
						int start = starts.get(k);
						if( ++i[k] < subSeqs.get(start).size() )
							break;
						i[k] = 0;
					}
					if( k < 0 )
						break; // break while
				}
				
			} // end of else
		} // end of for-loop
		
		return seqs;
	}
	
	/**
	 * (Multi-threaded)
	 * NOTE: Severe Efficiency drop! ???
	 * 
	 * Generate all of the possible sequence given the input nucleotide sequence
	 * and restrictions.
	 * 
	 * @param letters
	 *            A nucleotide sequence
	 * @param restrictions
	 *            Used to restrict the possible output sequence
	 * @return A collection of all of the RNA Sequences that satisfy the given
	 *         restrictions over the given nucleotide sequence. The elements in
	 *         the collection will not be in any particular order, but iterating
	 *         over the sequence will always happen in the same order.
	 */
	public static Collection<RNASequence> possibleSequencesMT(Letter[] letters, Restrictions restrictions) {
		int length = letters.length;
		
		int minMutationCount = Math.max( 0, restrictions.minMutationCount );
		int maxMutationCount = Math.min( length, restrictions.maxMutationCount );
		int minPairCount = Math.max( 0, restrictions.minPairCount );
		int maxPairCount = Math.min( length / 2, restrictions.maxPairCount );
		Set<Integer> validPairLocations = restrictions.validPairLocations;
		Set<Integer> validMutationLocations = restrictions.validMutationLocations;
		
		if( validPairLocations == null ) {
			validPairLocations = new HashSet<>();
			for( int i = 0; i < length; i++ )
				validPairLocations.add(i);
		}
		if( validMutationLocations == null ) {
			validMutationLocations = new HashSet<>();
			for( int i = 0; i < length; i++ )
				validMutationLocations.add(i);
		}
		
		// calculate all the possible mutation sets of the given length
		ArrayList<Integer> possibleMutationLocations = new ArrayList<>(length);
		for( int i = 0; i < length; i++ )
			if( validMutationLocations.contains(i) )
				possibleMutationLocations.add(i);
		Collection<Set<Integer>> mutationSets = subsets(possibleMutationLocations, minMutationCount, maxMutationCount);
		
		LinkedList<RNASequence> seqs = new LinkedList<>();
		ExecutorService threadPool = Executors.newFixedThreadPool( mutationSets.size() );
		for( Set<Integer> mutationSet : mutationSets ) {
			final Set<Integer> validPairLocationsF = validPairLocations;
			threadPool.submit( () -> {
				if( mutationSet.size() == 0 ) {
					// BASE CASE (for recursion): don't split the sequence

					// calculate the set of all possible pairing locations
					ArrayList<Pair> possiblePairs = new ArrayList<>(length * length - length);
					for( int n1 = 0; n1 < length; n1++ )
						for( int n2 = n1 + 1; n2 < length; n2++ )
							if( validPairLocationsF.contains(n1) && validPairLocationsF.contains(n2) && letters[n1].isCompatible(letters[n2]) ) {
								possiblePairs.add( new Pair(n1, n2) );
								if( restrictions.useOrderedPairs )
									possiblePairs.add( new Pair(n2, n1) );
							}

					// calculate all the valid subsets of possible pairs (ignoring mutations)
					Collection<Set<Pair>> pairSets = subsets(possiblePairs, minPairCount, maxPairCount, (set) -> {
						boolean[] flags = new boolean[length];
						for( Pair pair : set ) {
							if( flags[pair.getStart()] || flags[pair.getEnd()] )
								return false;
							flags[pair.getStart()] = true;
							flags[pair.getEnd()] = true;
						}
						return true;
					});

					// add each pair,mutation set combination to the collection of valid sequences
					for( Set<Pair> pairSet : pairSets ) {
						RNASequence rna = new RNASequence(letters, mutationSet, pairSet);
						synchronized(seqs) {
							seqs.add(rna);
						}
					}

				} else {  // mutationSet.size() != 0
					// split the sequence at each mutation since binds can not cross over mutations

					ArrayList<Integer> ends = new ArrayList<>( mutationSet.size() );
					ends.addAll(mutationSet);
					ends.add(length);

					TreeMap<Integer, ArrayList<RNASequence>> subSeqs = new TreeMap<>();

					{	int start = 0;
					for( int end : ends ) {
						// grab a sub-sequence between start (inclusive) and end (exclusive), and recurse
						Letter[] subLetter = new Letter[end - start];
						System.arraycopy(letters, start, subLetter, 0, subLetter.length);

						Restrictions subRestrictions = new Restrictions();
						subRestrictions.minMutationCount = 0;
						subRestrictions.maxMutationCount = 0;
						subRestrictions.minPairCount = 0;
						subRestrictions.maxPairCount = maxPairCount;
						subRestrictions.validMutationLocations = new HashSet<>();
						subRestrictions.validPairLocations  = new HashSet<>();
						for( int loc : validPairLocationsF ) {
							int newLoc = loc - start;
							if( 0 <= newLoc && newLoc < subLetter.length )
								subRestrictions.validPairLocations.add(newLoc);
						}
						subRestrictions.useOrderedPairs = restrictions.useOrderedPairs;

						Collection<RNASequence> c = possibleSequencesMT(subLetter, subRestrictions);
						subSeqs.put( start, new ArrayList<>(c) );

						// iterate
						start = end + 1;
					}
					}

					// for all valid combinations of 'subSeqs', join them into one RNASequence and add to the collection of valid sequence, 'seqs'
					ArrayList<Integer> starts = new ArrayList<>( subSeqs.keySet() );
					int[] i = new int[starts.size()];
					while(true) {
						// create RNASequence if valid, and add it to 'seqs'; otherwise, skip
						Set<Pair> pairSet = new HashSet<>();
						for( int k = 0; k < i.length; k++ ) {
							int start = starts.get(k);
							Set<Pair> subPairs = subSeqs.get(start).get(i[k]).getPairs();
							for( Pair p : subPairs )
								pairSet.add( p.shift(start) );
						}

						if( minPairCount <= pairSet.size() && pairSet.size() <= maxPairCount ) {
							RNASequence rna = new RNASequence(letters, mutationSet, pairSet);
							synchronized (seqs) {
								seqs.add(rna);	
							}
						}

						// iterate
						int k;
						for( k = i.length - 1; k >= 0; k-- ) {
							int start = starts.get(k);
							if( ++i[k] < subSeqs.get(start).size() )
								break;
							i[k] = 0;
						}
						if( k < 0 )
							break; // break while
					}

				} // end of else
			}); // end of lambda function (Runnable)
		} // end of for-loop
		try {
			threadPool.shutdown();
			while( !threadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS) ) {}
		} catch (InterruptedException ex) {
			throw new Error(ex); // for debugging
			// return seqs; // should not happen!
		}
		
		return seqs;
	}
	

	public int getLength() {
		return nodes.length;
	}
	
	public Node getNode(int position) {
		return nodes[position];
	}
	
	public Set<Pair> getPairs() {
		return pairs;
	}
	
	public boolean equals(Object obj) {
		if( this == obj )
			return true;
		else if( !(obj instanceof RNASequence) )
			return false;
		RNASequence seq = (RNASequence) obj;
		for( int i = 0; i < getLength(); i++ )
			if( !getNode(i).equals(seq.getNode(i)) )
				return false;
		return pairs.equals(seq.pairs);
	}
	
	public int hashCode() {
		return Arrays.hashCode(nodes) + 2147483647 * pairs.hashCode();
	}
	
	public String toString() {
		StringBuilder str = new StringBuilder();
		for( Node node : nodes )
			str.append(node);
		str.append(" ").append(pairs);
		return str.toString();
	}
	
}
