"""Parallel chunked HTTP range download with curl (resumable)."""
import subprocess, os, sys, concurrent.futures as cf, time
url, out, size = sys.argv[1], sys.argv[2], int(sys.argv[3])
CH = 32 * 1024 * 1024
tmp = out + ".chunks"; os.makedirs(tmp, exist_ok=True)


def get(i):
    a = i * CH; b = min(size, a + CH) - 1; p = os.path.join(tmp, f"{i:05d}")
    for attempt in range(30):
        if os.path.exists(p) and os.path.getsize(p) == b - a + 1:
            return i
        subprocess.run(["curl.exe", "--doh-url", "https://1.1.1.1/dns-query", "-sS", "--max-time", "900",
                        "-r", f"{a}-{b}", "-o", p, url], capture_output=True)
        time.sleep(2)
    raise RuntimeError(f"chunk {i} failed")


n = (size + CH - 1) // CH; t = time.time(); done = 0
with cf.ThreadPoolExecutor(16) as ex:
    for i in ex.map(get, range(n)):
        done += 1
        if done % 10 == 0:
            print(f"{done}/{n} {round(time.time() - t)}s", flush=True)
with open(out + ".part", "wb") as f:
    for i in range(n):
        with open(os.path.join(tmp, f"{i:05d}"), "rb") as c:
            f.write(c.read())
assert os.path.getsize(out + ".part") == size
os.replace(out + ".part", out)
for x in os.listdir(tmp):
    os.remove(os.path.join(tmp, x))
os.rmdir(tmp)
print("OK", out, flush=True)
