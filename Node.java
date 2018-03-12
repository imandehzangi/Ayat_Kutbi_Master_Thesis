package rna_sequences;

// immutable
public final class Node {
	
	private RNASequence seq;
	private int pos;
	private Letter letter;
	private boolean isMutated;
	private Pair pair;
	
	Node(RNASequence seq, int pos, Letter letter, boolean isMutated) {
		this.seq = seq;
		this.pos = pos;
		this.letter = letter;
		this.isMutated = isMutated;
		this.pair = null;
	}
	
	void setPair(Pair pair) {
		this.pair = pair;
	}
	
	public int getPosition() {
		return pos;
	}
	
	public Letter getLetter() {
		return letter;
	}
	
	public boolean isMutated() {
		return isMutated;
	}
	
	public Pair getPair() {
		return pair;
	}
	
	public boolean isStart() {
		if( pair == null )
			return false;
		return pair.getStart() == this.pos;
	}
	
	public boolean isEnd() {
		if( pair == null )
			return false;
		return pair.getEnd() == this.pos;
	}
	
	public boolean equals(Object obj) {
		if( this == obj )
			return true;
		if( !(obj instanceof Node) )
			return false;
		Node n = (Node) obj;
		return seq == n.seq
		    && pos == n.pos
		    && letter == n.letter
		    && isMutated == n.isMutated
		    && pair.equals(n.pair);
	}
	
	public int hashCode() {
		final int P = 61;
		
		int h = System.identityHashCode(seq);
		h = h * P + pos;
		h = h * P + letter.hashCode();
		h = h * P + pair.hashCode();
		if( isMutated )
			h = h * P + 1;
		
		return h;
	}
	
	public String toString() {
		return isMutated ? letter.toString().toLowerCase() : letter.toString().toUpperCase();
	}
	
}