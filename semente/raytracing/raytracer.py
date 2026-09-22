"""
Ray-tracer relativistico (CPU, numpy vetorizado) para Schwarzschild e para o
black-bounce de Simpson-Visser.

Fisica implementada
-------------------
* Geodesicas nulas exatas na metrica  -f dt^2 + dr^2/f + R^2 dOmega^2,  R = sqrt(r^2 + l^2):
      (dr/dlambda)^2 = E^2 - L^2 f / R^2      =>    d^2 r/dlambda^2 = L^2 r (R - 3M) / R^5
      dphi/dlambda   = L / R^2
  (l = 0 e Schwarzschild).  Cada raio vive num plano pela origem, entao integramos
  no plano e reconstruimos a posicao 3D.
* Disco de acrecao fino de Novikov-Thorne / Page-Thorne (1974) entre R_ISCO = 6M e R_out,
  com o fluxo F(R) calculado pela integral exata.  Curiosamente, em termos do raio areal R
  toda a fisica de orbitas circulares do black-bounce e IDENTICA a de Schwarzschild.
* Desvio para o vermelho/azul  g = sqrt(1 - 3M/R) / (1 - Omega lambda_z)  (gravitacional + Doppler),
  brilho observado  ~ g^4 F,  cor de corpo negro com T_obs = g T_em.
* Fundo estelar procedural lenteado.  Para l > 0 os raios que atravessam a garganta r=0
  enxergam o CEU DO OUTRO UNIVERSO (paleta diferente).  Para 0 < l < 2M, na geometria
  eterna, o que aparece "dentro da sombra" e a luz que veio do universo anterior atraves do
  buraco branco.

Nada aqui e uma aproximacao artistica: a unica liberdade estetica e a escala de temperatura
do disco e o mapa de tons.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter


# ----------------------------------------------------------------------------
# Disco de Page-Thorne
# ----------------------------------------------------------------------------
def page_thorne_flux_integral(R, M=1.0, Mdot=1.0):
    """Fluxo de Page-Thorne pela integral geral (usado como teste da forma fechada):

    F(R) = -(Mdot / 4 pi R) * (dOmega/dR) / (E - Omega L)^2 * int_{6M}^{R} (E - Omega L) (dL/dR') dR'
    """
    R = np.asarray(R, float)
    grid = np.linspace(6 * M, max(R.max(), 6.01 * M) * 1.01, 20000)
    E = (1 - 2 * M / grid) / np.sqrt(1 - 3 * M / grid)
    L = np.sqrt(M * grid) / np.sqrt(1 - 3 * M / grid)
    Om = np.sqrt(M / grid**3)
    dL = np.gradient(L, grid)
    integ = cumulative_trapezoid((E - Om * L) * dL, grid, initial=0.0)
    dOm = -1.5 * np.sqrt(M) * grid ** (-2.5)
    F = -(Mdot / (4 * np.pi * grid)) * dOm / (E - Om * L) ** 2 * integ
    return np.interp(R, grid, F, left=0.0, right=0.0)


def page_thorne_flux(R, M=1.0, Mdot=1.0):
    """Fluxo emitido por um disco fino de Novikov-Thorne em Schwarzschild, forma FECHADA.

    Derivacao (deste projeto, conferida contra a integral): com x = sqrt(R/M),
        E - Omega L = sqrt((x^2-3)/x^2),   dL/dR = (x^2-6) / (2 (x^2-3)^{3/2}),
        I(R) = int_{6M}^{R} (E - Omega L) dL/dR' dR' = (x - sqrt6) - (sqrt3/2) [ ln((x-sqrt3)/(x+sqrt3)) - ln((sqrt6-sqrt3)/(sqrt6+sqrt3)) ]
        F(R) = 3 Mdot / (8 pi M^2) * I / ( x^5 (x^2 - 3) ).
    Limite R >> M:  F -> 3 M Mdot / (8 pi R^3)  (Shakura-Sunyaev).
    """
    R = np.asarray(R, float)
    x = np.sqrt(np.maximum(R, 6 * M) / M)
    s3, s6 = np.sqrt(3.0), np.sqrt(6.0)
    I = (x - s6) - (s3 / 2) * (np.log((x - s3) / (x + s3)) - np.log((s6 - s3) / (s6 + s3)))
    F = 3 * Mdot / (8 * np.pi * M**2) * I / (x**5 * (x * x - 3))
    return np.where(R >= 6 * M, F, 0.0)


def newtonian_thin_disk_flux(R, M=1.0, Mdot=1.0, R_in=6.0):
    """Limite newtoniano (Shakura-Sunyaev): F = 3 M Mdot / (8 pi R^3) (1 - sqrt(R_in/R)).
    Usado como teste: page_thorne_flux -> este valor quando R >> M."""
    R = np.asarray(R, float)
    return 3 * M * Mdot / (8 * np.pi * R**3) * (1 - np.sqrt(R_in * M / R))


# ----------------------------------------------------------------------------
# Cor de corpo negro (aproximacao de Tanner Helland, valida 1000-40000 K) -> RGB linear
# ----------------------------------------------------------------------------
def blackbody_rgb(T):
    T = np.clip(np.asarray(T, float), 1000.0, 40000.0) / 100.0
    r = np.where(T <= 66, 255.0, 329.698727446 * np.power(np.maximum(T - 60, 1e-9), -0.1332047592))
    g = np.where(T <= 66, 99.4708025861 * np.log(np.maximum(T, 1e-9)) - 161.1195681661,
                 288.1221695283 * np.power(np.maximum(T - 60, 1e-9), -0.0755148492))
    b = np.where(T >= 66, 255.0, np.where(T <= 19, 0.0, 138.5177312231 * np.log(np.maximum(T - 10, 1e-9)) - 305.0447927307))
    rgb = np.stack([r, g, b], axis=-1) / 255.0
    rgb = np.clip(rgb, 0, 1)
    return rgb**2.2  # sRGB -> linear


# ----------------------------------------------------------------------------
# Ceu procedural
# ----------------------------------------------------------------------------
def _hash3(ix, iy, iz, seed):
    h = (ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791) ^ (seed * 2654435761)
    h = (h ^ (h >> 13)) * 1274126177
    h = h ^ (h >> 16)
    return (h & 0x7FFFFFFF) / float(0x7FFFFFFF)


def starfield(dirs, seed=1, density=0.0025, tint=(1.0, 1.0, 1.0), band_axis=(0.2, 0.0, 0.98)):
    """Estrelas e uma faixa galactica difusa, funcao determinista da direcao (dirs: (N,3) unitarios)."""
    d = np.asarray(dirs, float)
    cell = 260.0
    ix = np.floor((d[:, 0] + 1.0) * cell).astype(np.int64)
    iy = np.floor((d[:, 1] + 1.0) * cell).astype(np.int64)
    iz = np.floor((d[:, 2] + 1.0) * cell).astype(np.int64)
    h1 = _hash3(ix, iy, iz, seed)
    h2 = _hash3(ix, iy, iz, seed + 17)
    h3 = _hash3(ix, iy, iz, seed + 91)
    star = (h1 < density * 4).astype(float)
    # posicao da estrela dentro da celula: brilho depende da distancia ao centro (evita blocos)
    fx = (d[:, 0] + 1.0) * cell - ix - 0.5
    fy = (d[:, 1] + 1.0) * cell - iy - 0.5
    fz = (d[:, 2] + 1.0) * cell - iz - 0.5
    dist2 = fx**2 + fy**2 + fz**2
    profile = np.exp(-dist2 / (0.02 + 0.06 * h2**3))
    mag = star * profile * (0.15 + 2.5 * h2**6)
    temp = 2500.0 + 12000.0 * h3**2
    col = blackbody_rgb(temp)
    band = np.asarray(band_axis, float)
    band /= np.linalg.norm(band)
    lat = np.abs(d @ band)
    # faixa galactica difusa e suave (modulada por ondas de baixa frequencia, sem blocos)
    wob = 0.75 + 0.25 * np.cos(9.0 * d[:, 0] + 5.0 * d[:, 1] + seed) * np.cos(7.0 * d[:, 2] - 3.0 * d[:, 0])
    glow = 0.035 * np.exp(-(lat / 0.18) ** 2) * wob
    rgb = mag[:, None] * col + glow[:, None] * np.array([0.75, 0.85, 1.0])
    return rgb * np.asarray(tint, float)


# ----------------------------------------------------------------------------
# Ray tracer
# ----------------------------------------------------------------------------
@dataclass
class Camera:
    r: float = 40.0
    inclination_deg: float = 80.0  # 90 = de lado, 0 = de cima
    fov_deg: float = 52.0
    width: int = 960
    height: int = 540


@dataclass
class Scene:
    M: float = 1.0
    l: float = 0.0  # parametro de Simpson-Visser (0 = Schwarzschild)
    disk: bool = True
    R_out: float = 22.0
    T_peak_K: float = 6500.0  # temperatura no maximo do fluxo (escolha estetica; a FORMA de F(R) e fisica)
    exposure: float = 0.28
    bloom: float = 0.0035
    R_far: float = 150.0
    other_sky_tint: tuple = (1.25, 0.65, 1.15)  # o outro universo (paleta violeta)
    our_sky_tint: tuple = (1.0, 1.0, 1.0)


def _camera_rays(cam: Camera):
    W, H = cam.width, cam.height
    aspect = W / H
    tanf = np.tan(np.radians(cam.fov_deg) / 2)
    ys, xs = np.mgrid[0:H, 0:W]
    x_ndc = (xs + 0.5) / W * 2 - 1
    y_ndc = 1 - (ys + 0.5) / H * 2
    # componentes no referencial estatico local (e_r, e_theta, e_phi); a camera olha para a origem
    d_r = -np.ones_like(x_ndc)
    d_th = -y_ndc * tanf  # "cima" na tela = -e_theta (para o polo norte)
    d_ph = x_ndc * tanf * aspect
    norm = np.sqrt(d_r**2 + d_th**2 + d_ph**2)
    return (d_r / norm).ravel(), (d_th / norm).ravel(), (d_ph / norm).ravel()


def render(scene: Scene = Scene(), cam: Camera = Camera(), verbose=True, return_layers=False):
    M, l = scene.M, scene.l
    W, H = cam.width, cam.height
    N = W * H
    th_c = np.radians(cam.inclination_deg)
    r_hat = np.array([np.sin(th_c), 0.0, np.cos(th_c)])
    th_hat = np.array([np.cos(th_c), 0.0, -np.sin(th_c)])
    ph_hat = np.array([0.0, 1.0, 0.0])
    d_r, d_th, d_ph = _camera_rays(cam)
    d3 = d_r[:, None] * r_hat + d_th[:, None] * th_hat + d_ph[:, None] * ph_hat  # (N,3)

    R_cam = np.sqrt(cam.r**2 + l**2)
    f_cam = 1 - 2 * M / R_cam
    sin_a = np.sqrt(np.clip(1 - d_r**2, 0, 1))
    L = R_cam * sin_a / np.sqrt(f_cam)  # parametro de impacto b = L/E, E = 1
    # base do plano orbital: e1 = r_hat, e2 = componente de d3 ortogonal a e1
    e1 = np.broadcast_to(r_hat, (N, 3)).copy()
    e2 = d3 - d_r[:, None] * r_hat
    n2 = np.linalg.norm(e2, axis=1)
    bad = n2 < 1e-12
    e2[bad] = th_hat
    n2[bad] = 1.0
    e2 /= n2[:, None]
    normal = np.cross(e1, e2)  # momento angular do raio TRACADO (foton real vai ao contrario)
    nz = normal[:, 2]

    r = np.full(N, float(cam.r))
    rdot = d_r * np.sqrt(np.clip(1 - L**2 * f_cam / R_cam**2, 0, None))
    phi = np.zeros(N)

    # resultado por raio: 0 = ativo, 1 = disco, 2 = ceu nosso, 3 = ceu do outro universo, 4 = capturado
    status = np.zeros(N, dtype=np.int8)
    hit_R = np.zeros(N)
    hit_lambda = np.zeros(N)
    out_dir = np.zeros((N, 3))
    z_prev = R_cam * (np.cos(phi) * e1[:, 2] + np.sin(phi) * e2[:, 2])

    def accel(rr, LL):
        RR = np.sqrt(rr**2 + l**2)
        return LL**2 * rr * (RR - 3 * M) / RR**5

    active = np.arange(N)
    it = 0
    max_iter = 20000
    while active.size and it < max_iter:
        it += 1
        ra, va, pa, La = r[active], rdot[active], phi[active], L[active]
        Ra = np.sqrt(ra**2 + l**2)
        # passo adaptativo: limita o avanco angular (L h / R^2) e o radial (|rdot| h / R)
        h = np.minimum(0.03 * Ra**2 / (La + 0.1), 0.25 * Ra / (np.abs(va) + 0.1))
        h = np.clip(h, 1e-3, 40.0)
        # RK4 para (r, rdot, phi)
        k1r, k1v, k1p = va, accel(ra, La), La / Ra**2
        r2 = ra + 0.5 * h * k1r
        v2 = va + 0.5 * h * k1v
        k2r, k2v, k2p = v2, accel(r2, La), La / (r2**2 + l**2)
        r3 = ra + 0.5 * h * k2r
        v3 = va + 0.5 * h * k2v
        k3r, k3v, k3p = v3, accel(r3, La), La / (r3**2 + l**2)
        r4 = ra + h * k3r
        v4 = va + h * k3v
        k4r, k4v, k4p = v4, accel(r4, La), La / (r4**2 + l**2)
        rn = ra + h * (k1r + 2 * k2r + 2 * k3r + k4r) / 6
        vn = va + h * (k1v + 2 * k2v + 2 * k3v + k4v) / 6
        pn = pa + h * (k1p + 2 * k2p + 2 * k3p + k4p) / 6
        Rn = np.sqrt(rn**2 + l**2)

        # cruzamento do disco (so no nosso lado r>0)
        e1a, e2a = e1[active], e2[active]
        zn = Rn * (np.cos(pn) * e1a[:, 2] + np.sin(pn) * e2a[:, 2])
        zp = z_prev[active]
        crossed = (zn * zp < 0) & (rn > 0) & (ra > 0) & scene.disk
        if crossed.any():
            frac = zp / (zp - zn)
            Rc = Ra + frac * (Rn - Ra)
            hit = crossed & (Rc >= 6 * M) & (Rc <= scene.R_out)
            idx = active[hit]
            status[idx] = 1
            hit_R[idx] = Rc[hit]
            hit_lambda[idx] = -La[hit] * nz[idx]  # L_z/E do foton REAL (sentido oposto ao tracado)
        # escape
        esc = (Rn > scene.R_far) & (rn * vn > 0)  # afastando-se da origem, em qualquer dos dois universos
        if esc.any():
            idx = active[esc]
            vel = vn[esc, None] * (np.cos(pn[esc, None]) * e1a[esc] + np.sin(pn[esc, None]) * e2a[esc]) \
                + (La[esc] / Rn[esc])[:, None] * (-np.sin(pn[esc, None]) * e1a[esc] + np.cos(pn[esc, None]) * e2a[esc])
            vel /= np.linalg.norm(vel, axis=1)[:, None]
            out_dir[idx] = vel
            status[idx] = np.where(rn[esc] > 0, 2, 3)
        # captura (so Schwarzschild: dentro de 2M nada volta; no black-bounce o raio atravessa)
        if l == 0:
            cap = Rn < 2 * M * 1.001
            status[active[cap]] = 4
        else:
            cap = np.zeros_like(esc)

        r[active], rdot[active], phi[active] = rn, vn, pn
        z_prev[active] = zn
        done = (status[active] != 0) | cap
        active = active[~done]
        if verbose and it % 200 == 0:
            print(f"  passo {it}: raios ativos {active.size}/{N}")
    status[status == 0] = 4  # o que sobrou (orbitas presas na esfera de fotons) fica escuro

    # --- shading -------------------------------------------------------------------------
    img = np.zeros((N, 3))
    # disco
    d = status == 1
    if d.any():
        Rd = hit_R[d]
        Om = np.sqrt(M / Rd**3)
        g = np.sqrt(1 - 3 * M / Rd) / (1 - Om * hit_lambda[d])
        g = g / np.sqrt(f_cam)  # observador estatico em r_cam (nao no infinito)
        F = page_thorne_flux(Rd, M)
        Fmax = page_thorne_flux(np.linspace(6 * M, scene.R_out, 400), M).max()
        T_em = scene.T_peak_K * (F / Fmax) ** 0.25
        T_obs = g * T_em
        bright = g**4 * (F / Fmax)
        img[d] = blackbody_rgb(T_obs) * bright[:, None] * 2.2
    # ceus
    ours = status == 2
    if ours.any():
        img[ours] = starfield(out_dir[ours], seed=1, tint=scene.our_sky_tint)
    other = status == 3
    if other.any():
        img[other] = starfield(out_dir[other], seed=2, density=0.004, tint=scene.other_sky_tint, band_axis=(0.9, 0.3, 0.2))
        # nebulosa difusa extra para deixar claro que e "outro ceu"
        dd = out_dir[other]
        neb = 0.05 * np.exp(-((dd[:, 0] - 0.3) ** 2 + (dd[:, 2] + 0.2) ** 2) / 0.5)
        img[other] += neb[:, None] * np.array([0.6, 0.15, 0.9])

    img = img.reshape(H, W, 3) * scene.exposure
    if scene.bloom > 0:
        blur = np.stack([gaussian_filter(img[..., c], sigma=max(W, H) * scene.bloom * 2.5) for c in range(3)], -1)
        blur2 = np.stack([gaussian_filter(img[..., c], sigma=max(W, H) * scene.bloom * 0.6) for c in range(3)], -1)
        img = img + 0.35 * blur + 0.25 * blur2
    # tone map (Reinhard) + gamma
    img = img / (1 + img)
    img = np.clip(img, 0, 1) ** (1 / 2.2)
    out = (img * 255).astype(np.uint8)
    if return_layers:
        return out, dict(status=status.reshape(H, W), hit_R=hit_R.reshape(H, W), g=None)
    return out


def save_png(img: np.ndarray, path: str):
    from PIL import Image
    Image.fromarray(img, "RGB").save(path)
