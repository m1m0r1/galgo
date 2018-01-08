#ifndef LOGSUMEXP_H
#define LOGSUMEXP_H

float flogsumexp(const float* __restrict__ buf, int N);
double dlogsumexp(const double* __restrict__ buf, int N);

#endif // LOGSUMEXP_H
