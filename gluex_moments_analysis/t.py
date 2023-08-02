import gluex_moments_analysis.analyzeMomentsMatrix as ana
import numpy as np
import random

def get_Kmatrix():
  global Kmatrix
  try:
    return Kmatrix
  except:
    pass
  Kmatrix = ana.make_Kmatrix()
  return Kmatrix

def get_Kpseudoinverse():
  global Kinverse
  try:
    return Kinverse
  except:
    pass
  K = get_Kmatrix()
  u,e,vt = np.linalg.svd(K)
  Kinverse = np.zeros([K.shape[1], K.shape[0]], dtype=float)
  for i in range(len(e)):
    Kinverse[i,i] = 1 / e[i]
  Kinverse = vt.T @ Kinverse @ u.T

  # make a list of SU(N) generators out of the eigenvectors of K
  global Ksigma
  Ksigma = []
  for i in range(K.shape[1]):
    Ksigma.insert(0, ana.get_model1_rmatrix(rvector=vt[i,:]) / 2**0.5)
    rv = vt[i,:]
    rm = ana.get_model1_rmatrix(rvector=rv)
    rdiff = np.linalg.norm(rv - ana.get_model1_rvector(rmatrix=rm))
    if rdiff > 1e-12:
      print("bad rdiff", rdiff)
      return
  return Kinverse

def get_Bmatrix():
  global Bmatrix
  try:
    return Bmatrix
  except:
    pass
  Bmatrix = ana.make_Bmatrix()
  return Bmatrix

def get_SUNgenerators(from_Kmatrix=1):
  global SUNgenerators
  get_Kpseudoinverse()
  if from_Kmatrix:
    return Ksigma
  try:
    return SUNgenerators
  except:
    pass
  N = Ksigma[0].shape[0]
  SUNgenerators = []
  global rvectors
  rvectors = np.empty([N**2-1, N**2], dtype=float)
  # start off with N-1 diagonal generators
  for n in range(1, N):
    sigma = [1 for i in range(n)]
    sigma.append(-n)
    sigma += [0 for i in range(n+1, N)]
    sigma = np.array(sigma, dtype=complex)
    sigma /= np.linalg.norm(sigma) * 2**0.5
    sigma = np.diag(sigma)
    SUNgenerators.append(sigma)
    rvectors[n-1] = ana.get_model1_rvector(rmatrix=sigma)
  # next use Gramm-Schmidt to find a random set for the rest
  for n in range(N-1, N**2-1):
    r = np.random.uniform(-1, 1, [N**2])
    r[N**2 - 1] = -sum([r[i*N+i] for i in range(N-1)])
    for i in range(n):
      r = r - 2 * (r @ rvectors[i]) * rvectors[i]
    r /= np.linalg.norm(r) * 2**0.5
    SUNgenerators.append(ana.get_model1_rmatrix(rvector=r))
    rvectors[n] = r
  # test for orthogonality of the final rvectors
  unity = rvectors @ rvectors.T
  for i in range(N**2-1):
    if abs(unity[i,i] - 0.5) > 1e-10:
      print(f"disunity issue: unity({i},{i})={unity[i,i]}")
      return
    for j in range(i+1, N**2-1):
      if abs(unity[i,j]) > 1e-10:
        print(f"disunity concern: unity({i},{j})={unity[i,j]}")
        return
  SUNgenerators.append(np.diag(np.ones([N], dtype=complex) / (2*N)**0.5))
  return SUNgenerators

def test_Kmatrix(nrandom=1000):
  for i in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    h = ana.get_model1_moments(mX, abst)
    r = ana.get_model1_rvector(mX, abst)
    diff = np.linalg.norm(h - K @ r)
    print(f"{i}: {diff}")
    if diff > 1e-12:
      input("ok?")

def test_rmatrix(nrandom=1000):
  for i in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    rvector = ana.get_model1_rvector(mX, abst)
    rho = ana.get_model1_rmatrix(mX, abst)
    rrho = ana.get_model1_rmatrix(rvector=rvector)
    diff = np.linalg.norm(rrho - rrho)
    rbr = rvector @ b @ rvector
    tracediff = np.linalg.norm(np.trace(rho)**2 - np.trace(rho @ rho))
    print(f"{i}: {diff} {rbr} {tracediff}")
    if diff > 1e-12 or abs(rbr) > 1e-10 or tracediff > 1e-10:
      input("ok?")

def test_null_rmatrix():
  u,e,vt = np.linalg.svd(K)
  for i in range(e.shape[0], vt.shape[0]):
    h = K @ vt[i,:]
    hnorm = np.linalg.norm(h)
    print(f"{i}: {hnorm}")
    if hnorm > 1e-12:
      input("ok?")

def rank_of_rmatrix(nrandom=1000):
  Kinv = get_Kpseudoinverse()
  for i in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    rho1 = ana.get_model1_rmatrix(mX, abst)
    h = ana.get_model1_moments(mX, abst)
    rfromh = Kinv @ h
    rho2 = ana.get_model1_rmatrix(rvector=rfromh)
    u1,e1,vt1 = np.linalg.svd(rho1)
    print(f"{i}: ", e1[:3])
    u2,e2,vt2 = np.linalg.svd(rho2)
    rho1sqr = rho1 @ rho1
    rho2sqr = rho2 @ rho2
    print("   ", e2[:3], np.linalg.norm(rho1 - rho2),
                         np.linalg.norm(rho1sqr * np.trace(rho1sqr) - rho1sqr @ rho1sqr),
                         np.linalg.norm(rho2sqr * np.trace(rho2sqr) - rho2sqr @ rho2sqr))

def dig_for_rmatrix(nrandom=1000, niter=100, nalpha=207, truestart=0, match_check=0):
  Kinv = get_Kpseudoinverse()
  sigma = get_SUNgenerators()
  for itry in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    h = ana.get_model1_moments(mX, abst)
    rvector0 = Kinv @ h
    rho0 = ana.get_model1_rmatrix(rvector=rvector0)
    tracerho0 = np.real(np.trace(rho0))
    rho0sqr = rho0 @ rho0
    truerho = ana.get_model1_rmatrix(mX, abst)
    alpha = np.zeros(nalpha, dtype=float)
    a = np.zeros([len(sigma), nalpha], dtype=float)
    b = np.zeros([len(sigma)], dtype=float)
    w = np.zeros([nalpha, nalpha], dtype=float)
    if truestart > 0:
      for i in range(nalpha):
        alpha[i] = (np.real(np.trace(sigma[i] @ truerho)) -
                    np.real(np.trace(sigma[i] @ rho0)))
    for it in range(niter):
      rho = rho0
      for i in range(nalpha):
        rho = rho + alpha[i] * sigma[i]
      rhosqr = rho @ rho
      print(f"{it}: |alpha|=", np.linalg.norm(alpha), 
            ", |R rho0 - rho0**2|=", np.linalg.norm(rho0 * tracerho0 - rho0sqr),
            ", |R rho0 - rho**2|=", np.linalg.norm(rho * tracerho0 - rhosqr))
      for i in range(a.shape[0]):
        for j in range(a.shape[1]):
          tracerho_ij = np.trace(sigma[i] @ sigma[j] @ rho)
          tracerho0_ij = np.trace(sigma[i] @ sigma[j] @ rho0)
          a[i,j] = -np.real(tracerho_ij + tracerho0_ij)
          if i < nalpha:
            w[i,j] = 2 * np.real(tracerho_ij) / tracerho0
        if i < nalpha:
          a[i,i] = a[i,i] + tracerho0
          w[i,i] = w[i,i] - 1
        b[i] = (np.real(np.trace(sigma[i] @ rho0sqr)) -
                np.real(np.trace(sigma[i] @ rho0)) * tracerho0)
        if match_check:
          ci = [np.trace(sigma[i] @ sigma[j] @ rho0) for j in range(nalpha)]
          di = [np.trace(sigma[i] @ rho @ sigma[j]) for j in range(nalpha)]
          if i < a.shape[1]:
            print(f"{i} match check:", np.linalg.norm(tracerho0 * alpha[i]
                                       - ci @ alpha - di @ alpha - b[i]))
          else:
            print(f"{i} match check:", np.linalg.norm(0
                                       - ci @ alpha - di @ alpha - b[i]))
      ua,ea,vat = np.linalg.svd(a)
      eainv = np.zeros([a.shape[1], a.shape[0]], dtype=complex)
      eainv[:nalpha,:nalpha] = np.diag(1 / ea)
      alnew = vat.T @ eainv @ ua.T @ b
      alphadiff = alnew - alpha
      dalpha = alphadiff / 2
      print("   diff,step=", np.linalg.norm(alphadiff), np.linalg.norm(dalpha))
      ans = input("q to quit:")
      if ans == 'q':
        break
      alpha = alpha + dalpha
    rfinal = rvector0
    for i in range(nalpha):
      rfinal = rfinal + alpha[i] * vt[vt.shape[1] - i - 1,:]
    hfinal = K @ rfinal
    print("final result:", np.linalg.norm(rho * tracerho0 - rho @ rho), np.linalg.norm(hfinal - h))
    efinal,vtfinal = np.linalg.eigh(rho)
    print("final density spectrum is", efinal)
    """
    print("comparing a alpha with Tr(sigma rho0sqr)")
    a_alpha = a @ alpha
    for i in range(len(alpha)):
      print(f"{i}:", a_alpha[i], np.trace(sigma[i] @ rho0sqr))
    """

def seek_for_rmatrix(nrandom=1000, niter=100, truestart=0, eps=1):
  Kinv = get_Kpseudoinverse()
  sigma = get_SUNgenerators()
  for itry in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    h = ana.get_model1_moments(mX, abst)
    rvector0 = Kinv @ h
    rho0 = ana.get_model1_rmatrix(rvector=rvector0)
    tracerho0 = np.trace(rho0)
    rho0sqr = rho0 @ rho0
    truerho = ana.get_model1_rmatrix(mX, abst)
    alpha = np.zeros([len(sigma)], dtype=complex)
    a = np.zeros([len(sigma), len(sigma)], dtype=complex)
    b = np.zeros([len(sigma)], dtype=complex)
    if truestart:
      for i in range(len(sigma)):
        alpha[i] = np.trace(sigma[i] @ truerho)
    for it in range(niter):
      rho = rho0
      for i in range(len(sigma)):
        rho = rho + alpha[i] * sigma[i]
      rhosqr = rho @ rho
      print(f"{it}: nonzero, nonzero, zero=", np.linalg.norm(alpha), np.linalg.norm(rho0 * tracerho0 - rho0sqr), np.linalg.norm(rho * tracerho0 - rhosqr))
      rhocub = rhosqr @ rho
      for i in range(len(sigma)):
        for j in range(i, len(sigma)):
          a[i,j] = (8 * np.trace(sigma[i] @ rho) * np.trace(sigma[j] @ rho)
                  - 8 * np.trace(sigma[i] @ sigma[j] @ rhosqr)
                  - 4 * np.trace(sigma[i] @ rho @ sigma[j] @ rho))
          a[j,i] = a[i,j]
        a[i,i] = a[i,i] + 4 * np.trace(rhosqr)
        b[i] = (4 * np.trace(sigma[i] @ rho) * np.trace(rhosqr)
              - 4 * np.trace(sigma[i] @ rhocub))
      print(f"{it}: gradient magnitude", np.linalg.norm(b), abs(tracerho0**2 * np.trace(rhosqr) - np.trace(rho @ rhocub)))
      ea,ua = np.linalg.eigh(a)
      if ea[0] < -1e-3:
        print("negative curvature in surface, first few eigenvalues:", ea[:5])
        ans = input("q to quit, <number> for new epsilon: ")
        if ans == 'q':
          break
        try:
          eps = float(ans)
        except:
          pass
      dalpha = np.linalg.inv(a) @ b
      print("   step size is", np.linalg.norm(dalpha))
      alpha = alpha - dalpha * eps
    rfinal = rvector0
    for i in range(len(sigma)):
      rfinal = rfinal + alpha[i] * vt[i + len(e),:]
    hfinal = K @ rfinal
    print("final result:", np.linalg.norm(rho * tracerho0 - rho @ rho), np.linalg.norm(hfinal - h))
    efinal,vtfinal = np.linalg.eigh(rho)
    print("final density spectrum is", efinal)

def fiddle(nrandom=1000, truestart=0):
  Kinv = get_Kpseudoinverse()
  sigma = get_SUNgenerators(from_Kmatrix=0)
  nsigma = len(sigma)
  nalpha = nsigma - 1
  for i in range(len(sigma)):
    nonzeros = 0
    for j in range(sigma[0].shape[0]):
      for k in range(j+1,sigma[0].shape[1]):
        if abs(sigma[i][j,k]) > 1e-10:
          nonzeros = 1
          break
      if nonzeros > 0:
        break
    if nonzeros == 0:
      print("diagonal sigma:", i)
  global d
  d = np.zeros([nalpha, nalpha, nalpha], dtype=float)
  global f
  f = np.zeros([nalpha, nalpha, nalpha], dtype=float)
  for i in range(nalpha):
    for j in range(i, nalpha):
      sigma_ij = sigma[i] @ sigma[j]
      for k in range(j, nalpha):
        tracesigma_ijk = np.trace(sigma_ij @ sigma[k])
        d[i,j,k] = d[j,k,i] = d[k,i,j] = \
        d[i,k,j] = d[k,j,i] = d[j,i,k] = 4 * np.real(tracesigma_ijk)
        f[i,j,k] = f[j,k,i] = f[k,i,j] = 4 * np.imag(tracesigma_ijk)
        f[i,k,j] = f[k,j,i] = f[j,i,k] = -f[i,j,k]
  ntotal = 0
  for i in range(nalpha):
    nzeros = -1
    for j in range(nalpha):
      nonzeros = 0
      for k in range(nalpha):
        if abs(f[i,j,k]) > 1e-10:
          nonzeros = 1
          break
      if nonzeros == 0:
        nzeros += 1
    if nzeros > 0:
      print(f"{i}: nzeros=", nzeros)
    ntotal += nzeros
  print("total commuting pairs:", ntotal)
  while True:
    ans = input("test of f,d relations, enter k1,k2:")
    try:
      k = [int(kstr) for kstr in ans.split(',')]
    except:
      break
    dd = df = ff = 0
    for i in range(nalpha):
      for j in range(nalpha):
        dd += d[i,j,k[0]] * d[i,j,k[1]]
        df += d[i,j,k[0]] * f[i,j,k[1]]
        ff += f[i,j,k[0]] * f[i,j,k[1]]
    print("dd,df,ff=", [round(this,6) for this in (dd,df,ff)])
  return d,f

  for itry in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    h = ana.get_model1_moments(mX, abst)
    rvector0 = Kinv @ h
    rho0 = ana.get_model1_rmatrix(rvector=rvector0)
    tracerho0 = np.trace(rho0)
    rho0sqr = rho0 @ rho0
    truerho = ana.get_model1_rmatrix(mX, abst)
    alpha = np.zeros([nsigma], dtype=float)
    if truestart > 0:
      for i in range(nalpha):
        alpha[i] = (np.real(np.trace(sigma[i] @ truerho)) -
                    np.real(np.trace(sigma[i] @ rho0)))
    break
    """
    for it in range(niter):
      rho = rho0
      for i in range(len(sigma)):
        rho = rho + alpha[i] * sigma[i]
      rhosqr = rho @ rho
      print(f"{it}: nonzero, nonzero, zero=", np.linalg.norm(alpha), np.linalg.norm(rho0 * tracerho0 - rho0sqr), np.linalg.norm(rho * tracerho0 - rhosqr))
      rhocub = rhosqr @ rho
      for i in range(len(sigma)):
        for j in range(i, len(sigma)):
          a[i,j] = (8 * np.trace(sigma[i] @ rho) * np.trace(sigma[j] @ rho)
                  - 8 * np.trace(sigma[i] @ sigma[j] @ rhosqr)
                  - 4 * np.trace(sigma[i] @ rho @ sigma[j] @ rho))
          a[j,i] = a[i,j]
        a[i,i] = a[i,i] + 4 * np.trace(rhosqr)
        b[i] = (4 * np.trace(sigma[i] @ rho) * np.trace(rhosqr)
              - 4 * np.trace(sigma[i] @ rhocub))
      print(f"{it}: gradient magnitude", np.linalg.norm(b), abs(tracerho0**2 * np.trace(rhosqr) - np.trace(rho @ rhocub)))
      ea,ua = np.linalg.eigh(a)
      if ea[0] < -1e-3:
        print("negative curvature in surface, first few eigenvalues:", ea[:5])
        ans = input("q to quit, <number> for new epsilon: ")
        if ans == 'q':
          break
        try:
          eps = float(ans)
        except:
          pass
      dalpha = np.linalg.inv(a) @ b
      print("   step size is", np.linalg.norm(dalpha))
      alpha = alpha - dalpha * eps
    rfinal = rvector0
    for i in range(len(sigma)):
      rfinal = rfinal + alpha[i] * vt[i + len(e),:]
    hfinal = K @ rfinal
    print("final result:", np.linalg.norm(rho * tracerho0 - rho @ rho), np.linalg.norm(hfinal - h))
    efinal,vtfinal = np.linalg.eigh(rho)
    print("final density spectrum is", efinal)
    """

def explore(nrandom=10, truestart=0):
  Kinv = get_Kpseudoinverse()
  sigma = get_SUNgenerators()
  nsigma = len(sigma)
  nalpha = nsigma - 1
  for itry in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    h = ana.get_model1_moments(mX, abst)
    rvector0 = Kinv @ h
    rho0 = ana.get_model1_rmatrix(rvector=rvector0)
    tracerho0 = np.trace(rho0)
    rho0sqr = rho0 @ rho0
    truerho = ana.get_model1_rmatrix(mX, abst)
    alpha = np.zeros([nsigma], dtype=float)
    normalpha = tracerho0 * (15/16)**0.5
    if truestart > 0:
      for i in range(nalpha):
        alpha[i] = np.real(np.trace(sigma[i] @ truerho))
      print("norm of true alpha is", np.linalg.norm(alpha) / normalpha)
    zeros = 0
    for i in range(nalpha):
      if abs(alpha[i]) < 1e-10:
        zeros += 1
    print("non-zero components of true alpha:", nalpha - zeros)
    # fiddling around, find alpha for rho=diag(1,0,0...)
    rho = np.zeros(rho0.shape, dtype=complex)
    rho[itry,itry] = tracerho0
    for i in range(nalpha):
      alpha[i] = np.real(np.trace(sigma[i] @ rho))
    print("norm of test alpha is", np.linalg.norm(alpha) / normalpha)
    #print("alpha:", np.round(alpha, 6))
    zeros = 0
    for i in range(nalpha):
      if abs(alpha[i]) < 1e-10:
        zeros += 1
    print("non-zero components of alpha:", nalpha - zeros)
      
