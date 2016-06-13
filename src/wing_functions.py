import math

# c(y)
def f_chord(y, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2):
    if y >= 0 and y <= b_k/2:
        A_ = 2*(c_k - c_r)/b_k; B_ = c_r;
        return A_*y + B_
    elif y > b_k/2 and y <= b/2:
        A_ = (c_t - c_k)/(b/2 - b_k/2); B_ = c_k;
        return A_*(y - b_k/2) + B_
    
def f_chord_1(y, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2):
    A_ = 2*(c_k - c_r)/b_k; B_ = c_r;
    return A_*y + B_

def f_chord_2(y, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2):
    A_ = (c_t - c_k)/(b/2 - b_k/2); B_ = c_k;
    return A_*(y - b_k/2) + B_

# X_le(y)
def f_Xle(y, b_k, b, Lambda_le_1, Lambda_le_2):
    if y >= 0 and y <= b_k/2:
        A_ = math.tan(Lambda_le_1);
        return A_*y
    elif y > b_k/2 and y <= b/2:
        A_ = math.tan(Lambda_le_2);
        return (b_k/2)*math.tan(Lambda_le_1) + A_*(y - b_k/2)
    
def f_Xle_1(y, b_k, b, Lambda_le_1, Lambda_le_2):
    A_ = math.tan(Lambda_le_1);
    return A_*y
    
def f_Xle_2(y, b_k, b, Lambda_le_1, Lambda_le_2):
    A_ = math.tan(Lambda_le_2);
    return (b_k/2)*math.tan(Lambda_le_1) + A_*(y - b_k/2)

# X_le(y) * c(y)
def f_Xle_c(y, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2):
    return f_Xle(y, b_k, b, Lambda_le_1, Lambda_le_2)*f_chord(y, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2)

# y * c(y)
def f_y_c(y, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2):
    return y*f_chord(y, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2)

# eps_g(y)
def f_twist(y, eps_k, eps_t, b_k, b):
    if y >= 0 and y <= b_k/2:
        A_ = 2*eps_k/b_k; B_ = 0;
        return A_*y + B_
    elif y > b_k/2 and y <= b/2:
        A_ = (eps_t - eps_k)/(b/2 - b_k/2); B_ = eps_k;
        return A_*(y - b_k/2) + B_
    
def f_twist_1(y, eps_k, eps_t, b_k, b):
    A_ = 2*eps_k/b_k; B_ = 0;
    return A_*y + B_

def f_twist_2(y, eps_k, eps_t, b_k, b):
    A_ = (eps_t - eps_k)/(b/2 - b_k/2); B_ = eps_k;
    return A_*(y - b_k/2) + B_

# alpha0l(y)
def f_alpha0l(y, alpha0l_r, alpha0l_k, alpha0l_t, b_k, b):
    if y >= 0 and y <= b_k/2:
        A_ = 2*(alpha0l_k - alpha0l_r)/b_k; B_ = alpha0l_r;
        return A_*y + B_
    elif y > b_k/2 and y <= b/2:
        A_ = (alpha0l_t - alpha0l_k)/(b/2 - b_k/2); B_ = alpha0l_k;
        return A_*(y - b_k/2) + B_
    
def f_alpha0l_1(y, alpha0l_r, alpha0l_k, alpha0l_t, b_k, b):
    A_ = 2*(alpha0l_k - alpha0l_r)/b_k; B_ = alpha0l_r;
    return A_*y + B_

def f_alpha0l_2(y, alpha0l_r, alpha0l_k, alpha0l_t, b_k, b):
    A_ = (alpha0l_t - alpha0l_k)/(b/2 - b_k/2); B_ = alpha0l_k;
    return A_*(y - b_k/2) + B_

# [ alpha_0l(y) - eps_g(y) ] * c(y)
def f_alpha0l_epsg_c(y, c_r, c_k, c_t, 
                     eps_k, eps_t,
                     alpha0l_r, alpha0l_k, alpha0l_t, 
                     b_k, b, Lambda_le_1, Lambda_le_2):
    return (f_alpha0l(y, alpha0l_r, alpha0l_k, alpha0l_t, b_k, b)
               - f_twist(y, eps_k, eps_t, b_k, b)
           )*f_chord(y, c_r, c_k, c_t, b_k, b, Lambda_le_1, Lambda_le_2)