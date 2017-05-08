HTMLStart(filename = "dynamic",echo = TRUE)
sqrt(pi)
x = rnorm(10)
x^2

myfunction = function(x) return(summary(x))
### try to fix the function
myfunction(x)
plot(x)
HTMLplot()


HTMLStop()
