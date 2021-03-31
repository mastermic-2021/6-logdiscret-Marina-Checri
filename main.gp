g = Mod(6, 682492462409094395392022581537473179285250139967739310024802121913471471);
A = 245036439927702828116237663546936021015004354074422410966568949608523157;

n = 682492462409094395392022581537473179285250139967739310024802121913471471;


/* ****************************************************************************** */
/* ******************************* Explications : ******************************* */
/* ****************************************************************************** */

/* On cherche à calculer le logarithme discret de A modulo g.
 * Pour ce faire, nous allons utiliser
 * la méthode Baby-step Giant-step présentée en cours.
 
 * Cependant, pour réussir à résoudre efficacement le problème du log discret,
 * nous allons d'abord décomposer le problème
 * à l'aide de la décomposition de Pohlig-Hellman.
 * il s'agit donc de "casser" le problème en traitant des instances de log discrets
 * dans des groupes de cardinal *premier*.
 
 * On utilisera ensuite l'algorithme Baby-step Giant-step, comme annoncé ci-dessus.
*/



/* ****************************************************************************** */
/* ******************************** Fonctions : ********************************* */
/* ****************************************************************************** */

/*
* Baby-step Giant-step :
* Choisir B=ceil(sqrt(n)) // partie entière inférieure de racine de n
* calculer la liste L des g^{a_0}​ pour 0≤a_0<B
* poser G=g^{-B}
* pour 0≤a_1<B0, dès que AG^{a_1} est dans L et vaut g^{a_0}, rendre a = a_0 + a_1*B
*/

/* Nota Bene : Comme recommandé par Marc lors d'une discussion,
 * nous préfererons utiliser une map plutôt qu'une liste.
 * Les fonctions pour les maps sont mapdelete, mapget, mapisdefined et mapput */

baby_step_giant_step(A, g, n) = {
            my(B, bb, bs, gs);
	    if(A == 1, return(0) );
	    B  = ceil(sqrt(n));
	    bs = g; gs = 1;
	    bb = Map();

	    \\baby steps
 	    mapput(bb, 1, 0); \\ bb[0]=1
 	    mapput(bb, g, 1); \\ bb[1]=g

	    for(i=2, B-1,
	       bs=bs*g;
	       mapput(bb, bs , i)); \\ bb[i] = bs
 
	    \\giant steps
	    gs = (g^B)^(-1);
 	    for(i=0, B+1,
	       my(result);
	       if (mapisdefined(bb, A, &result)==1, return (result + B*i) );
  	       A = A*gs;
	    );
}




/* Lemme chinois ~ Décomposition de Pohlig Hellman :
* Rappel des idées

* G=<g> groupe cyclique d'ordre n=dm
  Pour tout A= g^a dans G, on a a = log_{g^d}(A^d)[m].

* Pour n=p^e, on calcule les log discrets dans le groupe <g^{p^{e-1}}> d'ordre p.
*/


\\Traitement des diviseurs de n-1, n=p^e.
traitement_diviseurs(g, p, e, a)={
	   my(s, g_, A);
	   s=0;
	   g_ = g^(p^(e-1));
	   for(i=0, e-1,
	      A = (a*g^(-s))^(p^(e-i-1));
	      my(tmp);
	      tmp = baby_step_giant_step(A, g_, p);
	      s = s + tmp*(p^i)
	   );
	   return (s);
}



Pohlig_Hellman(g, A, n)={
        my(fact, size, tab, p, q, e, d, g_, a_);
	fact = factor(n-1);
	size = matsize(fact)[1];
	tab = vector(size);
	for(i=1, size,
		p = fact[i, 1];
		e = fact[i, 2];
		q = p^e;
		d = (n-1)/q;
		g_ = g^d;
		a_ =  Mod(A,n)^d;
		tab[i] = Mod(traitement_diviseurs(g_, p, e, lift(a_)), q);
	);
	result = lift(Mod(lift(chinese(tab)),n));
	return (result);
}




A= Mod(A,n);
print(Pohlig_Hellman(g, A, n));



