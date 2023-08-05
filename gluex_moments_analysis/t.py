import gluex_moments_analysis.analyzeMomentsMatrix as ana
import numpy as np
import random
import scipy

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

def get_reference_rmatrix(N=16):
  rmatrix = np.zeros([N, N], dtype=complex)
  rmatrix[0,0] = 1
  return rmatrix

def get_SUNgenerators(kind=0, N=0):
  global SUNgenerators
  get_Kpseudoinverse()
  if kind == 0:
    return Ksigma
  try:
    return SUNgenerators
  except:
    pass

  if N == 0:
    N = Ksigma[0].shape[0]
  SUNgenerators = np.zeros([N**2, N, N], dtype=complex)
  if kind == 1:
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
      SUNgenerators[n-1] = sigma
      rvectors[n-1] = ana.get_model1_rvector(rmatrix=sigma)
    # next use Gramm-Schmidt to find a random set for the rest
    for n in range(N-1, N**2-1):
      r = np.random.uniform(-1, 1, [N**2])
      r[N**2 - 1] = -sum([r[i*N+i] for i in range(N-1)])
      for i in range(n):
        r = r - 2 * (r @ rvectors[i]) * rvectors[i]
      r /= np.linalg.norm(r) * 2**0.5
      SUNgenerators[n] = ana.get_model1_rmatrix(rvector=r)
      rvectors[n] = r

  else:
    # taken from wikipedia article, Structure Constants
    for n in range(2, N+1):
      gn = n**2 - 1
      dnorm = (2*n*(n - 1))**0.5
      SUNgenerators[gn-1][n-1][n-1] = (1 - n) / dnorm
      for m in range(1, n):
        anm = n**2 + 2*(m - n) - 1
        bnm = anm + 1
        SUNgenerators[gn-1][m-1][m-1] += 1 / dnorm
        SUNgenerators[anm-1][m-1][n-1] = 0.5
        SUNgenerators[anm-1][n-1][m-1] = 0.5
        SUNgenerators[bnm-1][m-1][n-1] = -0.5j
        SUNgenerators[bnm-1][n-1][m-1] = 0.5j

  # to duplicate the behavior of get_Kmatrix, append a unit matrix
  SUNgenerators[N**2 - 1] =  np.diag(np.ones([N], dtype=complex) / (2*N)**0.5)

  #print("checking orthonormality")
  for i in range(N**2 - 1):
    trace = np.trace(SUNgenerators[i] @ SUNgenerators[i])
    if abs(trace - 0.5) > 1e-10:
      print(f"{i}: bad normalization", trace)
      print(np.round(SUNgenerators[i], 6))
      return
    for j in range(i+1, N**2):
      trace = np.trace(SUNgenerators[i] @ SUNgenerators[j])
      if abs(trace) > 1e-10:
        print(f"{i},{j}: bad orthogonality", trace)
        return
  return SUNgenerators

def get_SUNtensors(sigma):
  global dijk
  global fijk
  nalpha = len(sigma) - 1
  dijk = np.zeros([nalpha, nalpha, nalpha], dtype=float)
  fijk = np.zeros([nalpha, nalpha, nalpha], dtype=float)
  for i in range(nalpha):
    for j in range(i, nalpha):
      sigma_ij = sigma[i] @ sigma[j]
      for k in range(j, nalpha):
        tracesigma_ijk = np.trace(sigma_ij @ sigma[k])
        dijk[i,j,k] = dijk[j,k,i] = dijk[k,i,j] = \
        dijk[i,k,j] = dijk[k,j,i] = dijk[j,i,k] = 4 * np.real(tracesigma_ijk)
        fijk[i,j,k] = fijk[j,k,i] = fijk[k,i,j] = 4 * np.imag(tracesigma_ijk)
        fijk[i,k,j] = fijk[k,j,i] = fijk[j,i,k] = -fijk[i,j,k]
  return dijk, fijk

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
        alpha[i] = 2 * (np.real(np.trace(sigma[i] @ truerho)) -
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
    tracerho0 = np.real(np.trace(rho0))
    rho0sqr = rho0 @ rho0
    truerho = ana.get_model1_rmatrix(mX, abst)
    alpha = np.zeros([len(sigma)], dtype=complex)
    a = np.zeros([len(sigma), len(sigma)], dtype=complex)
    b = np.zeros([len(sigma)], dtype=complex)
    if truestart:
      for i in range(len(sigma)):
        alpha[i] = 2 * np.trace(sigma[i] @ truerho)
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

def fiddle(nrandom=1000, truestart=0, kind=0, N=0):
  Kinv = get_Kpseudoinverse()
  sigma = get_SUNgenerators(kind=kind, N=N)
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
  d,f = get_SUNtensors(sigma)
  zeros = np.ones([nalpha, nalpha], dtype=int)
  for i in range(nalpha):
    for j in range(nalpha):
      for k in range(j, nalpha):
        if abs(d[i,j,k]) > 1e-10:
          zeros[j,k] = 0
          zeros[k,j] = 0
  print("total anticommuting pairs:", (np.sum(zeros) - nalpha)//2)
  minsetsize = sigma[0].shape[0]
  foundsets = []
  setstats = {'topsize': 0}
  def seekmutualsets(bag, level=0, i1=0):
    for i in range(i1, nalpha):
      #if level < 3:
      #  [print(" ", end="") for i in range(level)]
      #  print("starting new row", i)
      keep = 1
      for ib in bag:
        if zeros[i,ib] == 0:
          keep = 0
          break
      if keep:
        bag.append(i)
        if len(bag) >= minsetsize:
          print("found new set", bag)
          foundsets.append([baggie for baggie in bag])
        elif len(bag) > setstats['topsize']:
          print("found set of size", len(bag), bag)
          setstats['topsize'] = len(bag)
        seekmutualsets(bag, level+1, i+1)
        bag.pop()
  seekmutualsets([])
  print("found", len(foundsets), "anticommuting sets with >=",
        minsetsize, "members, top anticommuting set size",
        max([len(aset) for aset in foundsets]))

def explore(nrandom=1, kind=0, truestart=0, niter=1000000000, interactive=1):
  Kinv = get_Kpseudoinverse()
  sigma = get_SUNgenerators(kind=kind)
  nsigma = len(sigma)
  N = Ksigma[0].shape[0]
  for itry in range(nrandom):
    abst = random.uniform(0,2)
    mX = random.uniform(0.6,2.4)
    h = ana.get_model1_moments(mX, abst)
    rvector0 = Kinv @ h
    rho0 = ana.get_model1_rmatrix(rvector=rvector0)
    tracerho0 = np.real(np.trace(rho0))
    rho0sqr = rho0 @ rho0
    rho1 = get_reference_rmatrix(sigma[0].shape[0]) * tracerho0
    truerho = ana.get_model1_rmatrix(mX, abst)
    nalpha = nsigma - 1
    alpha = np.zeros([nalpha], dtype=float)
    agoal = np.zeros([nalpha], dtype=float)
    for i in range(nalpha):
      agoal[i] = 2 * np.real(np.trace(sigma[i] @ truerho))
    if truestart:
      for i in range(nalpha):
        alpha[i] = 2 * np.real(np.trace(sigma[i] @ truerho))
    else:
      for i in range(nalpha):
        alpha[i] = 2 * np.real(np.trace(sigma[i] @ rho1))
    try:
      global sigma_kind
      if kind == sigma_kind:
        d,f = dijk, fijk
    except:
      sigma_kind = kind
      d,f = get_SUNtensors(sigma)
    direction = -1
    lookback = 10000000
    redirects = [0] * lookback
    distances = [0] * lookback
    axes = [i for i in range(207, nalpha)]
    goals = [i for i in range(207, nalpha)]
    domegas = np.zeros([lookback, nalpha], dtype=float)
    domega = np.zeros([nalpha], dtype=float)
    for itry in range(niter):
      rho = (np.diag(np.ones([N], dtype=complex) * tracerho0 / N) +
             np.einsum('i,ijk', alpha, sigma[:nalpha]))
      dgoal = agoal[goals] - alpha[goals]
      distance = np.linalg.norm(dgoal)
      distances[itry % lookback] = distance
      normalpha = np.linalg.norm(alpha)
      normdomega = np.linalg.norm(domega)
      if normdomega < 1e-6:
        jacob = np.einsum('i,ijk', alpha, f[:,goals][:,:,axes])
        jacou,jacoe,jacovt = np.linalg.svd(jacob)
        jacoeinv = np.zeros(jacob.T.shape, dtype=float)
        for i in range(len(jacoe)):
          if jacoe[i] > 1e-4:
            jacoeinv[i,i] = 1/jacoe[i]
          elif jacoe[i] > 1e-8:
            jacoeinv[i,i] = 1/jacoe[i]**0.5
          elif i == 0:
            jacoeinv[i,i] = 1
          else:
            jacoeinv[i,i] = 0.01
        domega = np.zeros([nalpha], dtype=float)
        domega[axes] = np.real(jacovt.T @ jacoeinv @ jacou.T @ dgoal)
        domega *= direction
        normdomega = np.linalg.norm(domega)
      if normdomega > 1:
        dostep = np.real(domega) / normdomega
      else:
        dostep = np.real(domega)
      domegas[itry % lookback] = np.zeros([nalpha], dtype=float)
      domegas[itry % lookback] = dostep
      fdostep = np.einsum('ijk,k', f, dostep)
      def expm(m, prec=1e-20):
        result = np.diag(np.ones(m.shape[0], dtype=complex))
        mscale = np.linalg.norm(m)
        if mscale < prec:
          return result
        nsteps = int(mscale) + 1
        mfactor = np.array(m / nsteps, dtype=complex)
        for deg in range(1,1000):
          result = result + mfactor
          mfactor = mfactor @ m / nsteps / (deg + 1)
          if np.linalg.norm(mfactor) < mscale * prec:
            break
        if nsteps == 1:
          return result
        n = 1
        results = {}
        results[n] = np.array(result)
        while n * 2 <= nsteps:
          n *= 2
          result = result @ result
          results[n] = np.array(result)
        nsteps -= n
        while nsteps > 0:
          n //= 2
          if nsteps >= n:
            result = result @ results[n]
            nsteps -= n
        return result
      #dalpha = scipy.linalg.expm(fdostep_product) @ alpha - alpha
      exp_fdostep = expm(fdostep)
      dalpha = exp_fdostep @ alpha - alpha
      print("|rank1|,|alpha|,|dalpha|,|domega|,|dostep|,|alpha-truealpha|,dist=",
            f"{np.linalg.norm(tracerho0 * rho - rho @ rho):12.6e}",
            f"{np.real(normalpha) / tracerho0:18.15f}",
            f"{np.linalg.norm(dalpha):9.6f}",
            f"{np.linalg.norm(domega):9.6f}",
            f"{np.linalg.norm(dostep):9.6f}",
            f"{np.linalg.norm(alpha-agoal):9.6f}",
            f"{distance}")
      now = itry % lookback
      if interactive:
        ans = input("i to inspect, c to continue, q to quit:")
        if ans == 'c':
          if now > 1:
            redirects[now-1] = 1
          interactive = False
        elif ans == 'q':
          break
        elif ans == 'i':
          print("|rho - truerho| =", np.linalg.norm(rho - truerho))
          one = np.diag(np.ones(fdostep.shape[0], dtype=complex))
          print("|exp_fdostep * exp_fdostep.T - 1| =", np.linalg.norm(exp_fdostep @ exp_fdostep.T - one))
          print("exp_fdostep =\n", np.round(np.real(exp_fdostep), 6))
          print(len(axes), "axes:", axes)
          print(len(goals), "goals:", goals)
          print("jacou.T @ dgoal =\n", np.round(np.real(jacou.T @ dgoal), 6))
          print("jacoe=\n", np.round(jacoe, 6))
          print("jacovt @ alpha=\n", np.round(np.real(jacovt @ alpha[axes]), 6))
          while True:
            ans = input("dial:angle/+/-/s/q?")
            try:
              if len(ans) > 0 and ans[0] == '+':
                for axis in ans[1:].split(','):
                  if not int(axis) in axes:
                    axes.append(int(axis))
                print(len(axes), "axes:", axes)
                domega *= 0
                break
              elif len(ans) > 0 and ans[0] == '-':
                for axis in ans[1:].split(','):
                  if int(axis) in axes:
                    axes.remove(int(axis))
                print(len(axes), "axes:", axes)
                domega *= 0
                break
              sans = ans.split(':')
              dial = int(sans[0])
              angle = float(sans[1])
            except:
              if ans == 's':
                dalpha = trial_dalpha
                domega = trial_domega
              break
            trial_domega = np.zeros([nalpha], dtype=float)
            dials = [i for i in range(nalpha)]
            for axis in axes:
              dials.remove(axis)
              dials.insert(0,-axis)
            if dial < len(axes):
              trial_domega[axes] = np.real(jacovt[dial]) * angle
              print(f"twisting skewed axis {dial} by {angle}")
            else:
              trial_domega[dials[dial]] = angle
              print(f"twisting basis axis {dials[dial]} by {angle}")
            trial_fdomega = np.einsum('ijk,k', f, trial_domega)
            trial_exp_fdomega = expm(trial_fdomega)
            trial_alpha = trial_exp_fdomega @ alpha
            trial_dalpha = trial_alpha - alpha
            trial_rho = (np.diag(np.ones([N], dtype=complex) * tracerho0 / N) +
                         np.einsum('i,ijk', trial_alpha, sigma[:nalpha]))
            trial_dgoal = agoal[goals] - trial_alpha[goals]
            trial_distance = np.linalg.norm(trial_dgoal)
            trial_normalpha = np.linalg.norm(trial_alpha)
            print("|rank1|,|alpha|,|dalpha|,|domega|,|alpha-truealpha|,dist=",
                  f"{np.linalg.norm(tracerho0 * trial_rho - trial_rho @ trial_rho):12.6e}",
                  f"{np.real(trial_normalpha) / tracerho0:18.15f}",
                  f"{np.linalg.norm(trial_dalpha):9.6f}",
                  f"{np.linalg.norm(trial_domega):9.6f}",
                  f"{np.linalg.norm(trial_alpha - agoal):9.6f}",
                  f"{trial_distance}")
        domega -= dostep
      else:
        domega *= 0
      alpha = alpha + dalpha
      redirects[now] = 0
      if now > 5:
        if sum([redirects[i] for i in range(now-5, now)]) > 0:
          continue
        elif sum([distances[i+1] > distances[i] for i in range(now-5, now)]) == 5:
          print("******WARNING******* ", end="")
          print("regressing, try reversing direction")
          redirects[now] = 1
          direction *= -1
          interactive = 1
          continue
      if now > 10:
        if abs(distances[now] / distances[now-10] - 1) < 1e-5:
          if abs(distances[now]) < 1e-10:
            break
          print("******WARNING******* ", end="")
          print("stagnated, let's randomize the goals and go interactive")
          interactive = 1
          continue
        elif sum([np.linalg.norm(domegas[i+1] - domegas[i]) for i in range(now-10, now)]) < 1e-4:
          print("******WARNING******* ", end="")
          print("cycling, let's try a half-way point")
          interactive = 1
          continue
      
  rho = np.diag(np.ones([N], dtype=complex) * tracerho0/N)
  for i in range(nalpha):
    rho += alpha[i] * sigma[i]
  print("final difference: ", np.linalg.norm(rho - truerho))
  print("final rank-1 check: |R rho - rho @ rho| is ", np.linalg.norm(tracerho0 * rho - rho @ rho))
