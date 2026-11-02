import numpy as np
from PIL import Image
from collections import deque
im = np.asarray(Image.open("fig22-000.jpg").convert("RGB")).astype(int)[:, :470]
r, g, b = im[..., 0], im[..., 1], im[..., 2]
masks = {"S_blue": (b > 150) & (r < 90) & (g < 130),
         "C_red": (r > 170) & (g < 90) & (b < 90),
         "QS_black": (r < 60) & (g < 60) & (b < 60)}
def blobs(mask, minpix, maxpix):
    lab = np.zeros(mask.shape, int); n = 0; out = []
    for y0, x0 in zip(*np.nonzero(mask)):
        if lab[y0, x0]: continue
        n += 1; q = deque([(y0, x0)]); lab[y0, x0] = n; pts = []
        while q:
            y, x = q.popleft(); pts.append((y, x))
            for dy, dx in ((1,0),(-1,0),(0,1),(0,-1)):
                yy, xx = y + dy, x + dx
                if 0 <= yy < mask.shape[0] and 0 <= xx < mask.shape[1] and mask[yy, xx] and not lab[yy, xx]:
                    lab[yy, xx] = n; q.append((yy, xx))
        if minpix <= len(pts) <= maxpix:
            p = np.array(pts); out.append((p[:, 0].mean(), p[:, 1].mean(), len(pts), np.ptp(p[:, 0]), np.ptp(p[:, 1])))
    return out
res = {}
for k, m in masks.items():
    bl = sorted(blobs(m, 25, 400), key=lambda t: t[1]); res[k] = bl
    print(k, [(round(y), round(x), n, int(h), int(w)) for y, x, n, h, w in bl])
gray = np.asarray(Image.open("fig22-000.jpg").convert("L")).astype(int)
prof = (gray[:, 30:62] < 110).sum(1); rows = np.nonzero(prof)[0]
grp, cur = [], [rows[0]]
for rr in rows[1:]:
    if rr - cur[-1] <= 2: cur.append(rr)
    else: grp.append(cur); cur = [rr]
grp.append(cur); print("left label rows:", [(round(np.mean(x), 1), len(x)) for x in grp])
colp = (gray[505:530, :470] < 110).sum(0); cols = np.nonzero(colp)[0]
grp, cur = [], [cols[0]]
for cc in cols[1:]:
    if cc - cur[-1] <= 3: cur.append(cc)
    else: grp.append(cur); cur = [cc]
grp.append(cur); print("bottom label cols:", [(round(np.mean(x), 1), len(x)) for x in grp])
import pickle; pickle.dump(res, open("dig32.pkl", "wb"))
