def encode_num_to_length(num, length):
  return ('%0' + str(2*length) + 'x') % int(num)
