package rna_sequences;

// immutable
public final class Pair {

	private int start, end;

	Pair(int start, int end) {
		this.start = start;
		this.end = end;
	}

	public int getStart() {
		return start;
	}
	
	public int getEnd() {
		return end;
	}
	
	public Pair shift(int d) {
		return new Pair(start + d, end + d);
	}
	
	public String toString() {
		return "(" + start + "," + end + ")";
	}
	
	public boolean equals(Object obj) {
		if( this == obj )
			return true;
		if( !(obj instanceof Pair) )
			return false;
		Pair p = (Pair) obj;
		return start == p.start && end == p.end;
	}
	
	public int hashCode() {
		return start + 131 * end;
	}

}