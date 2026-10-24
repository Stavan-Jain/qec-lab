from balance_literals import lit
def expand(block,depth,x):
 cap=1<<depth
 for k in range(depth-1,-1,-1):
  g=1<<k;mask=((1<<(2*block*cap))-1)//((1<<(2*block*g))-1)*((1<<(block*g))-1);x=(x|(x<<(block*g)))&mask
 return x
def packlit(xs,width):
 if not any(xs):return '0'
 cap=1<<(max(1,len(xs))-1).bit_length();depth=cap.bit_length()-1
 small=1
 while small<max([x.bit_length()for x in xs]+[1]):small*=2
 small=min(small,width)
 packed=sum(x<<(small*i)for i,x in enumerate(xs));s=lit(packed);b=small
 while b<width:
  s=f'(PackedDigits.expandData {b} {depth} {s})';packed=expand(b,depth,packed);b*=2
 assert packed==sum(x<<(width*i)for i,x in enumerate(xs))
 return s
