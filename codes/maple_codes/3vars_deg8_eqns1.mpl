alias(t=RootOf(x^3-2)):
sol:=solve([-b_0_11*t+b_0_12*t^2-2*b_0_13+2*b_0_14*t, -b_1_11*t-2*b_1_13+2*b_1_14*t-b_2_11*t^2-b_3_8*t^2-b_4_7*t^2-b_5_6*t^2-4*t^2, -b_1_13*t^2-b_2_11*t-2*b_2_13-b_3_9*t^2-b_4_8*t^2-b_5_7*t^2-2*b_5_9*t, -2*b_3_11*t-4*b_3_13+4*b_3_14*t-2*b_4_11*t^2-2*b_6_8*t^2-b_7_7*t^2, 2*b_3_14-b_4_11*t+b_4_12*t^2+2*b_4_14*t+2*b_5_12+2*b_7_9-12,2*b_4_14-b_5_11*t+b_5_12*t^2, b_6_12*t^2-2*b_6_13-2*b_7_13*t, b_6_12*t-b_6_13*t^2-2*b_7_13+2*b_7_14*t-8*t^2, 2*b_7_14-b_8_11*t+b_8_12*t^2+2*b_8_14*t+2*b_9_12, 2*b_8_14-b_9_11*t+b_9_12*t^2+8, -2*b_11_13+2*b_11_14*t-8*t, -2*b_11_13*t^2+2*b_11_14+2*b_12_14*t-8, -b_11_13*t-b_11_14*t^2+4*b_12_14+4*t^2, -b_11_14*t+b_12_14*t^2+4*t]):

# parse to giac
# map(x->rhs(x),sol);
# map(x->lhs(x),sol);
