Comments are in the code prefixed by @note

Very good code, but I would suggest to comment using Doxygen the major classes and functions
Moreover, you have used
```cpp
return std::move(matrix_C);
```
This is not only uselessbut may even degrade preformance. Just return the object by value. The compiler will optimize and possibly move since a return value is an rvalue.
So just write
```cpp
return matrix_C;
```
I tell you more precisely why. How is the returned value of you matrix-matrix multiplication used if you return by value? If it is used in a context where it is possible to move 
then the compiler will move the value. If it is used to copy the value, then the compiler will copy the value. So, in any case, the compiler will do the best thing for you. So, just return the value by value.
With the std::move you are not gaining anything. Moreover if you have 
```cpp
Matrix C=A*B;
```
then the compiler may implement copy elision (return value optimization) just by constructing C with the value provided by the multiplication (if is returns a value). If you use std::move the compiler will implement a move, which is however less efficient than copy elision!. Again, return by value is better.

By the way, the compiler should have given you a warning about the std::move.
```bash
Matrix_implementation.hpp: In instantiation of ‘algebra::Matrix<T, StorageOrder> algebra::operator*(const Matrix<T, StorageOrder>&, const Matrix<T, StorageOrder>&) [with T = double; Ordering StorageOrder = algebra::Ordering::RowMajor]’:
main.cpp:61:24:   required from here
Matrix_implementation.hpp:659:34: warning: moving a local object in a return statement prevents copy elision [-Wpessimizing-move]
  659 |         return std::move(matrix_C);
  ```

**A negative note** Matrix file is not given, so I could not run the test easily.
