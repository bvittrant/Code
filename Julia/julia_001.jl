########################### Test Julia ########################################

# 10 october 2018
# Benjamin Vittrant

# https://docs.julialang.org/en/v1/manual/mathematical-operations/

###############################################################################

# Some maths and print

x = 5
y = 10
z = x + y
println(z)

# Some math constant

pi

###############################################################################

## characters
### (https://en.wikibooks.org/wiki/Introducing_Julia/Strings_and_characters)

x = "Hello world"
x
y = "Hello"
z = " world!"
y*z

# Interpolation
x = 10
"The value of x is $(x)"

# Test variables types
typeof(1)

##############################################################################

# Updating operators

x = 1
x += 3
x

###############################################################################

# Vectorized "dot" operators

[1,2,3] .^ 3
print([1,2,3] .^ 3)
print([1,2,3] .* 3)
print([1,2,3] .+ 3)

f(x,y) = 3x + 4y
A = [1.0, 2.0, 3.0]
B = [4.0, 5.0, 6.0]
print(f.(pi,A))
print(f.(A,B))
print(f.(A,B).*10)
print(f.(A,B) .*10 .+1)

###############################################################################

# Array and Plot test
# https://docs.julialang.org/en/v1/manual/arrays/

d = randn((8,8))

using Pkg
Pkg.add("Gadfly")
using Gadfly
plot(y=[1,2,3])

using Pkg
Pkg.add("PyPlot")
using Pyplot
x = -2pi:0.1:2pi;
plot(x, sin(x.^2)./x);
