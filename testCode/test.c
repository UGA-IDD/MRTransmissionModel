



int a; // a is an int
int *b; // b is a pointer to an int
int **c; // c is a pointer to a pointer to an int


int i = 10; //i is an int, it has allocated storage to store an int.
int *k; // k is an uninitialized pointer to an int.
//It does not store an int, but a pointer to one.
k = &i; // make k point to i. We take the address of i and store it in k
int j = *k; //here we dereference the k pointer to get at the int value it points
//to. As it points to i, *k will get the value 10 and store it in j


