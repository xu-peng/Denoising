from pyirt import irt

src_fp = open(file_path,'r')

# alternatively, pass in list of tuples in the format of [(user_id, item_id, ans_boolean)]
# ans_boolean is 0/1.


# (1)Run by default
item_param, user_param = irt(src_fp)

# (2)Supply bounds
item_param, user_param = irt(src_fp, theta_bnds = [-4,4], alpha_bnds=[0.1,3], beta_bnds = [-3,3])

# (3)Supply guess parameter
guessParamDict = {1:{'c':0.0}, 2:{'c':0.25}}
item_param, user_param = irt(src_fp, in_guess_param = guessParamDict)