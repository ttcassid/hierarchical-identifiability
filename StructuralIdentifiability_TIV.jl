# StructuralIdentifiability for TIV model

using StructuralIdentifiability

# Define ODE model
ode1 = @ODEmodel(
	# T'(t) = λ - β*V(t)*T(t) - δₜ *T] # state variable
	# I'(t) = β*V(t)*T(t) - δ*I(t) # state variable
	# V'(t) = p*I(t) - c*V(t) # state variable
	x1'(t) = λ - β*x3(t)*x1(t) - δₜ*x1(t), # state variable
	x2'(t) = β*x3(t)*x1(t) - δ*x2(t), # state variable
	x3'(t) = p*x2(t) - c*x3(t), # state variable
	y1(t) = x3(t), # outputs
	# y2(t) = x2(t), # outputs
	# y3(t) = x3(t) # outputs
)

# assess identifiability
assess_identifiability(ode1)
