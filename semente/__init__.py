"""
SEMENTE - Cosmologia de Buracos Negros
======================================

Pacote open source que resolve numericamente, aprende neuralmente e desenha
a cadeia de ideias:

    colapso gravitacional  ->  buraco negro  ->  (ricochete)  ->  buraco branco
                           ->  Big Bang de um novo universo

Modulos
-------
geometry    : metricas (Schwarzschild, Simpson-Visser black-bounce), Kruskal, Penrose, Flamm
geodesics   : integrador geral de geodesicas (Christoffel via sympy) + forma orbital u(phi)
raytracer   : imagem fisicamente realista de buraco negro / ponte para o outro universo
collapse    : colapso de Oppenheimer-Snyder: o interior da estrela E um universo FRW
interior    : interior de Schwarzschild como cosmologia de Kantowski-Sachs
bounce      : ricochete cosmologico: Einstein-Cartan (torcao/spin) e LQC; universo-filho
selection   : selecao natural cosmologica (Smolin) - dinamica populacional de universos
pinn        : rede neural informada por fisica que aprende a solucao do ricochete
figures     : todas as figuras do projeto

Unidades: G = c = 1 (geometricas). Quando hbar aparece, usamos unidades de Planck.
"""
__version__ = "0.1.0"
