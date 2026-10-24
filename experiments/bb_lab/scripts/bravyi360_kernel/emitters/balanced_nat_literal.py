"""Render kernel naturals with bounded literal parse size and a balanced shift/OR DAG."""
import re,sys

def literal(n,chunk_bits=8192):
 if n.bit_length()<=chunk_bits:return hex(n)
 chunks=(n.bit_length()+chunk_bits-1)//chunk_bits
 shift=(chunks//2)*chunk_bits
 return f'(Nat.lor {literal(n&((1<<shift)-1),chunk_bits)} (Nat.shiftLeft {literal(n>>shift,chunk_bits)} {shift}))'

def replace(s,chunk_bits=8192):
 return re.sub(r'0x[0-9a-fA-F]{2049,}',lambda m:literal(int(m[0],16),chunk_bits),s)
if __name__=='__main__':
 src,dst=sys.argv[1:3];s=open(src).read();open(dst,'w').write(replace(s))
