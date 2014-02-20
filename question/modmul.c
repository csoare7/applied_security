#include "modmul.h"
#include "string.h"

#define noBits 1024

void readDevRand(){
	int bytes = 64;
	char data[bytes];
	FILE *fp;
	fp=fopen("/dev/random", "r");
	fread(data, 1, bytes, fp);
	fclose(fp);
}


mp_limb_t getOmega(mpz_t N){
	mp_limb_t t = 1;
	mp_limb_t limb_N = mpz_getlimbn(N, 0);

	for (int i = 1; i <= mp_bits_per_limb-1; i++){
		t = (t*t*limb_N); 
	}
	//gmp_printf ("%Mu\n", -t);
	return -t;
}

void getRhoSquared( mpz_t *t, mpz_t N){
	mpz_set_ui(*t,1);
	int lN = mpz_size(N);

	for (int i=1; i <= 2*lN*mp_bits_per_limb; i++){
		mpz_add(*t, *t, *t);
		mpz_mod(*t, *t, N);
	}
}

void MontMul(mpz_t *r, mpz_t x, mpz_t y, mpz_t N){
	int lN = mpz_size(N);
	
	mpz_t temp1, temp2;
	mp_limb_t r0, yi, x0, omega, u;
	
	mpz_set_ui(*r, 0);
	omega = getOmega(N);
	mpz_init(temp1);mpz_init(temp2);

	for(int i=0; i<lN; i++){

		r0 = mpz_getlimbn(*r, 0);
		//gmp_printf("r0= %Mu\n", r0);
		yi = mpz_getlimbn(y, i);
		x0 = mpz_getlimbn(x, 0);

		u = (r0 + yi * x0) * omega;
		
	    //r = (r + yi * x + u * N)/b;
		mpz_mul_ui(temp1, x, yi);
		mpz_mul_ui(temp2, N, u);
		mpz_add(temp1, temp1, temp2);
		mpz_add(*r, *r, temp1);

		int size = mpz_size(*r);
		mp_limb_t a[size];
		mpz_export(a, NULL, 1, sizeof(mp_limb_t), 1, 0, *r);
		mpz_import(*r, size-1, 1, sizeof(mp_limb_t), 1, 0, a);

	}
	if(mpz_cmp(*r, N) >= 0 ){
		mpz_sub(*r, *r, N);
	}
}

void MontgomeryMul(mpz_t *r, mpz_t x, mpz_t y, mpz_t N){
	mpz_t x_p, y_p, r_p, rho, one;
	mpz_init(x_p);mpz_init(y_p);mpz_init(rho);mpz_init(r_p);mpz_init(one);
	mpz_set_ui(one, 1);
	getRhoSquared(&rho, N);

	MontMul(&x_p, x, rho, N);
	MontMul(&y_p, y, rho, N);
	MontMul(&r_p, x_p, y_p, N);
	MontMul(r, r_p, one, N);
}

void MontExp(mpz_t *r, mpz_t x, mpz_t y, mpz_t N){
	mpz_t t_p, x_p, rho, one, temp;
	mpz_init(t_p); 
	mpz_init(x_p);
	mpz_init(rho);
	mpz_init(one);
	mpz_init(temp);

	getRhoSquared(&rho, N);
	mpz_set_ui(one, 1);
	
	//mpz_mod(temp, x, N);
	//mpz_set(x, temp);
	
	MontMul(&t_p, one, rho, N);
	MontMul(&x_p, x, rho, N);

	
	// gmp_printf("t_p %Zx\n",x_p);
	
	for(int i = noBits-1; i >= 0; i--){
		MontMul(r, t_p, t_p, N);
		//mpz_set(t_p, temp);
		if(mpz_tstbit(y, i)){
			MontMul(&t_p, *r, x_p, N);
		}
		else mpz_set(t_p, *r);
		//gmp_printf("i= %d\n t_p= %Zx\n",i , t_p);
	}
	MontMul(r, t_p, one, N);
}

int max(int a, int b){
	return ((a) > (b) ? (a) : (b));
}

void windowedExponent(mpz_t x, mpz_t y, double k, mpz_t N){
	int size = (int) pow(2, k)/2;
	mpz_t T[size]; 
	mpz_t temp, t;

	int l, u;

	//char s[1026];

	mpz_init(temp);
	mpz_init(t);
	
	////////////////////////////////////////////////////////////
	/////// Precompute T = [[i]x | i in {1, 3, .., 2^k - 1}] ///
	////////////////////////////////////////////////////////////

	mpz_init(T[0]);
  	mpz_set(T[0],x);

  for (int i=1; i<size; i++) {
    mpz_init(T[i]);
    mpz_mul(T[i],T[i-1],x);
    //MontgomeryMul(&(T[i]), T[i-1], x, N);
    mpz_mul(T[i],T[i],x);
    mpz_mod(T[i],T[i],N);
  }
	//mpz_clear(temp);

	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////

	mpz_set_ui(t, 1); // set t to neutral element of the group
	//mpz_get_str(s, 2, y); //convert y into base-2
	//printf("dan\n");
	//printf("%s\n", s);
	//printf("%d\n", i);
	//int i = 0;	
	// while (i <= (strlen(s)-1)){
	// 	//printf("%Zd\n", strlen(s)-1);
	// 	if (s[i] == '0'){
	// 		mpz_mul(t, t, t);
	// 		mpz_mod(t, t, N);
	// 		i++;
	// 	}
	// 	else{
	// 		l = min( (int) i+k-1, strlen(s)-1);
	// 		while(s[l] == '0'){
	// 			l--; 
	// 		}
	// 		for (int h = 0; h < l-i+1; h++){
	// 			mpz_mul(t, t, t);
	// 			mpz_mod(t, t, N);
	// 		}
	// 		int u = 0;
	// 		int e = 1;
	// 		for(int bit = l; bit >= i; bit--){
	// 			if(s[bit]=='1'){
	// 				u=u+e;
	// 			}
	// 			e=e<<1;
	// 		}
	// 		if (u != 0){
	// 			mpz_mul(t, t, T[(u-1/2)]);
	// 			mpz_mod(t, t, N);	
	// 		}
	// 		i = l + 1;
	// 	}
	// }
	// gmp_printf("%Zx\n", t);
	//mpz_clear(t);
	//mpz_clear(T);
	u = 0;
	int i = noBits - 1; // i = |y| - 1
	while (i >= 0){
		if(mpz_tstbit(y,i) == 0){
			l = i;
			u = 0;
			/*mpz_mul(t, t, t);
	 		mpz_mod(t, t, N);*/
		}
		else{
			l = max(i-k+1, 0);
			while(mpz_tstbit(y, l) == 0) l++;
			u = 0;
			int e = 1;
			for(int bit = l; bit <= i; bit++){
				if(mpz_tstbit(y, bit)){
					u=u+e;
				}
				e=e<<1;
			}
		}

		for (int h = 0; h < i-l+1; h++){
	 		//mpz_mul(t, t, t);
	 		//mpz_mod(t, t, N);
	 		MontgomeryMul(&t, t, t, N);
	 	}
    
		if (u != 0){
	 		mpz_mul(t, t, T[((u-1)/2)]);
	 		mpz_mod(t, t, N);	
	 	}
		i = l - 1;
	}
	gmp_printf("%Zx\n", t);
	mpz_clear(t);
	//mpz_clear(T);
}

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
		//mpz_powm(c, m, e, N);
		//windowedExponent(m, e, 6, N);
		MontExp(&c, m, e, N );

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

void stage5(){
	//mpz_t N;

	//mpz_init(N);

	//gmp_scanf("%Zd", N);


	//getOmega(N);
	//getRhoSquared(N);
	//getRhoSquared();
	//readDevRand();
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
   else if( !strcmp( argv[ 1 ], "stage5" ) ) {
    stage5();
  }

  return 0;
}

