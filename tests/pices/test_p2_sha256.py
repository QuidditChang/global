"""Verify the embedded SHA-256 against independent standard-library hashes."""
import hashlib,shutil,subprocess,tempfile,unittest
from pathlib import Path
class SHA256(unittest.TestCase):
    def test_vectors_and_chunk_boundaries(self):
        root=Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);c=d/'sha.c';c.write_text('#include <stdio.h>\n#include "pices_sha256.h"\nint main(void){PicesSHA s;char out[65];unsigned char b[37];size_t n;pices_sha_init(&s);while((n=fread(b,1,37,stdin)))pices_sha_add(&s,b,n);pices_sha_end(&s,out);puts(out);return 0;}\n')
            subprocess.run([shutil.which('cc'),'-I'+str(root/'lib'),str(c),'-o',str(d/'sha')],check=True)
            for data in [b'',b'abc',b'a'*55,b'a'*56,b'a'*63,b'a'*64,bytes(range(256))*4096]:
                got=subprocess.check_output([str(d/'sha')],input=data).decode().strip();self.assertEqual(got,hashlib.sha256(data).hexdigest())
if __name__=='__main__':unittest.main()
