# upa
Aquest és l'inici d'una aventura ...


Linear solvers: 

	Implemented: 

		- c++ LU decomp with partial pivoting  -> All non-singular matrixes work; overwrites initial matrix

		- c++ Cholesky decomposition           -> Symm, pos def matrixes; does preserve initial matrix

	Suggestions: 

		- f90 counterparts for LU, Chol

		- QR decomposition

		- Probably need 2 versions for each decomposition: 

		         One which overwrites the initial matrix (memory saving) + one that preserves it (debugging or other)
