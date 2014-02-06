#include "modmul.h"

/*
Perform stage 1:

- read each 3-tuple of N, e and m from stdin,
- compute the RSA encr yption c,
- then write the ciphertext c to stdout.
*/

void RSA_Encrypt(){
	mpz_t N, e, m, c;	
		
	//init
	mpz_init(N);
	mpz_init(e);
	mpz_init(m);
	mpz_init(c);
	
	while(gmp_scanf("%Zx", N) != EOF){
		//read
		gmp_scanf("%Zx", e);
		gmp_scanf("%Zx", m);
	
		//Encryption
		mpz_powm(c, m, e, N);

		gmp_printf("%Zx\n", c);
	}

	mpz_clear(N);
	mpz_clear(e);
	mpz_clear(m);
	mpz_clear(c);
		
}

void stage1() {
	RSA_Encrypt();
}

/*
Perform stage 2:

- read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
- compute the RSA decryption m,
- then write the plaintext m to stdout.
*/
void compute_CRT(mpz_t p, mpz_t q, mpz_t d_p, mpz_t d_q, mpz_t i_q, mpz_t c, mpz_t m){
	mpz_t m1, m2, h, temp1, temp2;

	mpz_init(m1);
	mpz_init(m2);
	mpz_init(h);
	mpz_init(temp1);
	mpz_init(temp2);

	mpz_powm(m1, c, d_p, p); //m1 = ( c^d_p ) ( mod p )
	mpz_powm(m2, c, d_q, q); //m2 = ( c^d_q ) ( mod q )

	mpz_sub(temp1, m1, m2); // temp1 = m1 - m2
	mpz_mul(temp2, i_q, temp1); // temp2 = i_q * (m1 - m2)
	mpz_mod(h, temp2, p); // temp2 mod p

	mpz_mul(temp1, h, q); //h*q

	mpz_add(m, m2, temp1); // m = m2 + hq

	mpz_clear(m1);
	mpz_clear(m2);
	mpz_clear(h);
	mpz_clear(temp1);
	mpz_clear(temp2);
}


void RSA_Decrypt(){
	mpz_t N, d, p, q, d_p, d_q, i_p, i_q, c, m ;

	mpz_init(N);
	mpz_init(d);
	mpz_init(p);
	mpz_init(q);
	mpz_init(d_p);
	mpz_init(d_q);
	mpz_init(i_p);
	mpz_init(i_q);
	mpz_init(c);
	mpz_init(m);
	
	while(gmp_scanf("%Zx", N) != EOF){
		//read
		gmp_scanf("%Zx", d);
		gmp_scanf("%Zx", p);
		gmp_scanf("%Zx", q);
		gmp_scanf("%Zx", d_p);
		gmp_scanf("%Zx", d_q);
		gmp_scanf("%Zx", i_p);
		gmp_scanf("%Zx", i_q);
		gmp_scanf("%Zx", c);
		//compute m
		//mpz_powm(m, c, d, N);

		compute_CRT(p, q, d_p, d_q, i_q, c, m);
		
		gmp_printf("%Zx\n", m);	
	}
	//clear
	mpz_clear(N);
	mpz_clear(d);
	mpz_clear(p);
	mpz_clear(q);
	mpz_clear(d_p);
	mpz_clear(d_q);
	mpz_clear(i_p);
	mpz_clear(i_q);
	mpz_clear(c);
	mpz_clear(m);
}

void stage2() {
	RSA_Decrypt();
}

/*
Perform stage 3:

- read each 5-tuple of p, q, g, h and m from stdin,
- compute the ElGamal encryption c = (c_1,c_2),
- then write the ciphertext c to stdout.
*/

void ElGamal_Encrypt(){
	mpz_t p, q, g, h, m, y, c1, c2, s, tmp;
	
	mpz_init(p); //large modulus
	mpz_init(q); //small modulus
	mpz_init(g); //generator of Fp with order q
	mpz_init(h); //public key
	mpz_init(m); //plaintext
	 
	mpz_init(y); //random key
	mpz_init(c1);
	mpz_init(c2);
	mpz_init(s);
	mpz_init(tmp);

	while(gmp_scanf("%Zx", p) != EOF){	
		gmp_scanf("%Zx", q);
		gmp_scanf("%Zx", g);
		gmp_scanf("%Zx", h);
		gmp_scanf("%Zx", m);

		//choose random y. Set y=1 for now
		mpz_set_str(y, "1",10);

		//compute c1
		mpz_set(c1, g);

		//compute c2
		mpz_set(s, h);
		mpz_mul(tmp, m, s);
		mpz_mod(c2, tmp, p);

		gmp_printf("%Zx\n", c1);
		gmp_printf("%Zx\n", c2);	
	}
	mpz_clear(p); //large modulus
	mpz_clear(q); //small modulus
	mpz_clear(g); //generator of Fp with order q
	mpz_clear(h); //public key
	mpz_clear(m); //plaintext
	 
	mpz_clear(y); //random key
	mpz_clear(c1);
	mpz_clear(c2);
	mpz_clear(s);
	mpz_clear(tmp);
}

void stage3() {
	ElGamal_Encrypt();
}

/*
Perform stage 4:

- read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
- compute the ElGamal decryption m,
- then write the plaintext m to stdout.
*/

void ElGamal_Decrypt(){
	mpz_t p, q, g, x, c1, c2, m, tmp, s;

	mpz_init(p); //large modulus
	mpz_init(q); //small modulus
	mpz_init(g); //generator of Fp with order q
	mpz_init(x); //private key
	mpz_init(c1); 
	mpz_init(c2); 
	mpz_init(m);
	mpz_init(tmp);
	mpz_init(s);

	while(gmp_scanf("%Zx", p) != EOF){	
	  	gmp_scanf("%Zx", q);
	  	gmp_scanf("%Zx", g);
	  	gmp_scanf("%Zx", x);
	  	gmp_scanf("%Zx", c1);
	  	gmp_scanf("%Zx", c2);

		mpz_powm(tmp, c1, x, p ); //compute s=c1^x
		mpz_invert(s, tmp, p); // invert s=s-1

		mpz_mul(tmp, c2, s);
		mpz_mod(m, tmp, p);

	  	gmp_printf("%Zx\n", m);	
	  }
	mpz_clear(p); //large modulus
	mpz_clear(q); //small modulus
	mpz_clear(g); //generator of Fp with order q
	mpz_clear(x); //private key
	mpz_clear(c1); 
	mpz_clear(c2); 
	mpz_clear(m);
	mpz_clear(tmp); 
	mpz_clear(s);
}

void stage4() {
	ElGamal_Decrypt();

}

/*
The main function acts as a driver for the assignment by simply invoking
the correct function for the requested stage.
*/

int main( int argc, char* argv[] ) {
  if     ( !strcmp( argv[ 1 ], "stage1" ) ) {
    stage1();
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    stage2();
  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    stage3();
  }
  else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    stage4();
  }

  return 0;
}
