package rna_sequences;

public enum Letter {
	
	A, C, G, U;
	
	public static Letter valueOf(char letter) {
		switch(letter) {
			case 'A': return A;
			case 'C': return C;
			case 'G': return G;
			case 'U': return U;
		}
		return null;
	}
	
	public static Letter[] seq(char[] letters) {
		Letter[] seq = new Letter[letters.length];
		for( int i = 0; i < seq.length; i++ )
			seq[i] = valueOf( letters[i] );
		return seq;
	}
	
	public static Letter[] seq(String letters) {
		return seq( letters.toCharArray() );
	}
	
	public boolean isCompatible(Letter letter) {
		return this.ordinal() + letter.ordinal() == 3;
	}
	
	public static String toString(Letter[] seq) {
		StringBuilder str = new StringBuilder();
		for( Letter x : seq )
			str.append( x.toString() );
		return str.toString();
	}
	
}
