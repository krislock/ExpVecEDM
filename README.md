## ExpVecEDM

ExpVecEDM is a MATLAB code for solving the Euclidean Distance Matrix completion
problem of finding locations of n points in r-dimensions that satisfy given
partial distance constraints.

`X = ExpVecEDM(Dpartial,A,r)` solves the Euclidean Distance Matrix completion
problem with partial distance (squared) matrix `Dpartial`, and anchor positions
`A`, in dimension `r`.

## Dependencies

- [CVX](http://cvxr.com/cvx/)

## Testing

In the ExpVecEDM directory, run the following command to test the code:

```
>> runtest
```

## License

The GNU Public License, Version 3.0+

## Citation Information

Dmitriy Drusvyatskiy, Nathan Krislock, Yuen-Lam Voronin, and Henry Wolkowicz.
Noisy Euclidean distance realization: Robust facial reduction and the pareto
frontier. SIAM Journal on Optimization, 27(4):2301-2331, 2017.
[https://doi.org/10.1137/15M103710X](https://doi.org/10.1137/15M103710X)


