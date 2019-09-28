# Test neural network in Julia

###############################################################################
# Library importation
using LinearAlgebra # For dot product

###############################################################################
# STEP 0
# Read input and output

# We create an input array
X = [1 0 1 0 ; 1 0 1 1 ; 0 1 0 1]
# We create an output array
Y = [1 ; 1 ; 0]

###############################################################################
# We create a sigmoide function what will be our activation function
# .(+-*/) syntax for matrix operation
sig(x) = 1 ./(1 .+ exp.(-x))
# ... And its derivative
d_sig(x) = x .*(1 .-x)

###############################################################################
# STEP 1
#  Initialize weights and biases with random values (There are methods to
# initialize weights and biases but for now initialize with random values)

# Initialization of variables
epoch = 5000 # Setting training iteration
lr = 0.1 # Learning rate
inputlayer_neuron = size(X)[2] # number of features in data set
hiddenlayer_neurons = 3 # number of features in data set
output_neurons = 1 # number of neurons at output layer

# Weight and bias initialization
wh = rand(inputlayer_neuron, hiddenlayer_neurons)
bh = rand(1, hiddenlayer_neurons)
wout = rand(hiddenlayer_neurons,output_neurons)
bout = rand(1, output_neurons)

###############################################################################
# Loop for neural network work
println(Y)
for i in collect(1:1:5000)
    global wh, bh, wout, bout

    ###########################################################################
    ### Forward propogation ###
    # STEP 2: Calculate hidden layer input
    hidden_layer_input1 = X * wh
    hidden_layer_input  = hidden_layer_input1 .+ bh
    # STEP 3: Perform non linear transformation on hidden linear input
    hiddenlayer_activations = sig(hidden_layer_input)
    # STEP 4: Step 4: Perform linear and non-linear transformation of hidden
    # layer activation at output laye
    output_layer_input1 = hiddenlayer_activations * wout
    output_layer_input = output_layer_input1 .+ bout
    output = sig(output_layer_input)
    println(output)

    ###########################################################################
    ### Backward propagation ###
    # STEP 5: Calculate gradient of Error(E) at output layer
    E = Y - output
    # STEP 6: Compute slope at output and hidden layer
    slope_output_layer = d_sig(output)
    slope_hidden_layer = d_sig(hiddenlayer_activations)
    # STEP 7: Compute delta at output layer
    d_output = E .* slope_output_layer
    # STEP 8: Calculate Error at hidden layer
    Error_at_hidden_layer = d_output * wout'
    # STEP 9: Compute delta at hidden layer
    d_hiddenlayer = Error_at_hidden_layer * slope_hidden_layer
    # STEP 10: Update weight at both output and hidden layer
    wout = wout .+ hiddenlayer_activations' * d_output .* lr
    bout = bout .+ sum(d_output, dims=(1)) .* lr
    # STEP 11: Update biases at both output and hidden layer
    wh = wh .+ X' * d_hiddenlayer .* lr
    bh = bh .+ sum(d_hiddenlayer, dims=(1)) .* lr
    end
