import subprocess

output = subprocess.run(["./NTS_WEBAPP_DD"], input="2\n1\n1.0\n0.95\n1\n1\n5\n1\n1\n5\n1\n1\n0\n0\n0\n0\n0.00001\n".encode())