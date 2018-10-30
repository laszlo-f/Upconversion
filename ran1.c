double ran1(){//random numbers on the interval [0,1]
  assert(2147483647==RAND_MAX); //warning:  only has precision 1/RAND_MAX
  return (double)(rand()) / (double)(RAND_MAX);
}
