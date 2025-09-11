# Generating matrices
##question1
A <- matrix(c(1,4,2,5,3,6), nrow = 3, byrow = TRUE)
B <- matrix(c(0,-2,-1,-3,1,-4), nrow = 3, byrow = TRUE)
A+B

##question2
A <- matrix(c(1,2,2,3), nrow = 2, byrow = TRUE)
B <- matrix(c(-1,0,1,1), nrow = 2, byrow = TRUE)
A-B

##question3
A <- matrix(c(1,2,0,3,-1,-2), nrow = 2, byrow = TRUE) 
C <- 3*A
C
##question4
A <- matrix(c(1,4,2,5,3,6), nrow = 3, byrow = TRUE)
B <- matrix(c(0,-2,-1,-3,1,-4), nrow = 3, byrow = TRUE)
C <- A %*% B 
C
##question5
A <- matrix(c(1,1,1), nrow = 3, byrow = TRUE)
B <- matrix(c(1,1,1), nrow = 1, byrow = TRUE)
C <- A %*% B 
C
# Note: %*% performs matrix multiplication (inner dimensions must match).
# The * operator is for element-wise multiplication (requires identical dimensions).
##question6 
A <- matrix(c(1,1,1), nrow = 1, byrow = TRUE)
B <- matrix(c(1,1,1), nrow = 3, byrow = TRUE)
C <- A %*% B
C
##question7 
A <- matrix(c(0,0,1,1,0,0,0,1,0), nrow = 3, byrow = TRUE)
B <- matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, byrow = TRUE)
C <- A %*% B 
C
##question8
A <- matrix(c(1,5,2,0,3,1), nrow = 3, byrow = TRUE)
B <- t(A)
B
#question9
A <- matrix(c(1,-1,0,2,3,2), nrow = 3, byrow = TRUE)
a_31 <- A[3, 1]
a_22 <- A[2, 2]
result <- a_31 - a_22
print(result)

#question10
# The diag() function is the easiest way to do the identity matrix
I <- diag(3)

# Note: R has a built-in constant for pi, and 3^-1 is written as 1/3.
A <- matrix(c(pi, 3, -1,
              2, 0, 0,
              5, 1/3, 0.5), nrow = 3, byrow = TRUE)
C <- I %*% A
print(C)

#question11
known_matrix <- matrix(c(1, 2, 
                         2, 3, 
                         0, 6), nrow = 3, byrow = TRUE)
y <- known_matrix[1, 2]
x <- known_matrix[3, 1]
print(paste("The value of x is:", x))
print(paste("The value of y is:", y))

#question12
A <- matrix(c(2, 1, 
              0, 1), nrow = 2, byrow = TRUE)
A_inverse <- solve(A)
print(A_inverse)

#question13
B <- matrix(c(3, 1, 3, 2), nrow = 2, byrow = TRUE)
determinant_B <- det(B) #Calculate the determinant using the det() function.
trace_B <- sum(diag(B)) #Calculate the trace by summing the diagonal elements.
print(paste("The determinant of B is:", determinant_B))
print(paste("The trace of B is:", trace_B))

#question16
A <-matrix(c(1, 1, 1, 1), nrow = 1, byrow = TRUE)
B <- matrix

#question19
A <- matrix(c(1,  1,  1,
              1, -1,  1,
              1, -1, -1), nrow = 3, byrow = TRUE)
b <- c(3, 3, -1)
x <- solve(A, b)#Solve the system for the vector 'x' using solve(A, b).
print(x)

#questions20
A <- matrix(c(1,  2, -1, -1,
              2, -1, -1,  1,
              2,  1,  3, -1,
              3,  3,  1,  1), nrow = 4, byrow = TRUE)

b <- c(-2, 6, -4, 1)
y <- solve(A, b) #solve system Ay = b for y using the solve() function.
print(y)








c <- 10*A[,1]
x <- B[2,2]

D <- cbind(A,c)
K <- rbind(A,c)

# Addition/subtraction

A+B
A-B
