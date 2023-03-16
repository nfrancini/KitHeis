import numpy as np
import pylab as pl
import sys
import os
from scipy.optimize import curve_fit

# FUNZIONI PER L'ANALISI DEGLI ERRORI
def jack_estimate_1obs(func, obs, block_size, index):
    resample = np.concatenate((obs[:index*block_size], obs[block_size*index+block_size:]))
    return func(resample)

def jack_error_1obs(func, obs, block_size):
    est = np.zeros(len(obs)//block_size)
    for x in range(0, len(obs)//block_size):
        est[x] = jack_estimate_1obs(func, obs, block_size, x)
    sum1 = np.sum(est**2)
    sum2 = np.sum(est)
    sum1 = sum1/(len(obs)/block_size)
    sum2 = (sum2/(len(obs)/block_size))**2
    return np.sqrt(len(obs)/block_size - 1) * np.sqrt(sum1-sum2)

def jack_estimate_2obs(func, obs1, obs2, block_size, index):
    resample1 = np.concatenate((obs1[:index*block_size], obs1[block_size*index+block_size:]))
    resample2 = np.concatenate((obs2[:index*block_size], obs2[block_size*index+block_size:]))
    return func(resample1, resample2)

def jack_error_2obs(func, obs1, obs2, block_size):
    est = np.zeros(len(obs1)//block_size)
    for x in range(0, len(obs1)//block_size):
        est[x] = jack_estimate_2obs(func, obs1, obs2, block_size, x)
    sum1 = np.sum(est**2)
    sum2 = np.sum(est)
    sum1 = sum1/(len(obs1)/block_size)
    sum2 = (sum2/(len(obs1)/block_size))**2
    return np.sqrt(len(obs1)/block_size - 1) * np.sqrt(sum1-sum2)

def jack_estimate_3obs(func, obs1, obs2, obs3, block_size, index):
    resample1 = np.concatenate((obs1[:index*block_size], obs1[block_size*index+block_size:]))
    resample2 = np.concatenate((obs2[:index*block_size], obs2[block_size*index+block_size:]))
    resample3 = np.concatenate((obs3[:index*block_size], obs3[block_size*index+block_size:]))
    return func(resample1, resample2, resample3)

def jack_error_3obs(func, obs1, obs2, obs3, block_size):
    est = np.zeros(len(obs1)//block_size)
    for x in range(0, len(obs1)//block_size):
        est[x] = jack_estimate_3obs(func, obs1, obs2, obs3, block_size, x)
    sum1 = np.sum(est**2)
    sum2 = np.sum(est)
    sum1 = sum1/(len(obs1)/block_size)
    sum2 = (sum2/(len(obs1)/block_size))**2
    return np.sqrt(len(obs1)/block_size - 1) * np.sqrt(sum1-sum2)

# FUNZIONI PER ESTRARRE LA GIUSTA OSSERVABILE
# MEDIA
def media(obs):
    return np.mean(obs)

# PARAMETRO DI BINDER (devo passare il modulo quadro del vettore in gioco)
def binder(obs):
    return 1-np.mean(obs**4)/(3*(np.mean(obs**2))**2)

# CALORE SPECIFICO (da capire le convenzioni su eventusle CT^2)
def spec_heat(obs):
    return (np.mean(obs**2) - np.mean(obs)**2)*V/T**2

# SUSCETTIVITÀ
def susc(obs_x, obs_y, obs_z):
    return ((np.mean(obs_x**2) - np.mean(obs_x)**2) + (np.mean(obs_y**2) - np.mean(obs_y)**2) + (np.mean(obs_z**2) - np.mean(obs_z)**2))*V/T
# def K2(obs):
#     return (np.mean(obs**2) - np.mean(obs)**2)*(V**2)
#
# def K3(obs):
#     return (np.mean(obs**3) - 3*np.mean(obs**2)*np.mean(obs) + 2*np.mean(obs)**3)*(V**3)
#
# def K4(obs):
#     return (np.mean(obs**4) - 4*np.mean(obs**3)*np.mean(obs) + 12*np.mean(obs**2)*(np.mean(obs)**2) - 3*np.mean(obs**2)**2 - 6*np.mean(obs)**4)*(V**4)
#
# def binder(obs):
#     return (np.mean(obs**2))/((np.mean(obs))**2)
#
# def corr_lenght(obs1, obs2):
#     return np.sqrt((np.mean(obs1)/np.mean(obs2)) -1) / (2*np.sin(np.pi/L))

# LETTURA DEL FILE DA TERMINALE
# PROVA A FARE BASEPATH DELLA CARTELLA E
# CICLARE SU TUTTI I FILE NELLA CARTELLA
basepath = sys.argv[1]
filenames = os.listdir(basepath)
path_filenames = []

for entry in filenames:
    path_filenames.append(os.path.join(basepath, entry))

for fname in path_filenames:
    print(fname)
    # LA PRIMA RIGA È FORMATA DA NUMERO DI CELLE, NUMERO DI SITI, NUMERO DI COMPONENTI DELLO SPIN, TEMPERATURA ED ALPHA
    L, V, K, T, alpha = np.genfromtxt(fname, dtype = "double", delimiter = "\t", unpack = True, max_rows = 1)
    ene_dens, mx, my, mz, nx, ny, nz, zx, zy, zz, sx, sy, sz = np.genfromtxt(fname, dtype = "double", delimiter = "\t", unpack = True, skip_header = 2, invalid_raise = False)

    # SCARTO LA TESTA DEL FILE SCARTO IL PRIMO DECIMO
    skip = len(ene_dens)//10

    ene_dens = ene_dens[skip:]
    mx = mx[skip:]
    my = my[skip:]
    mz = mz[skip:]
    nx = nx[skip:]
    ny = ny[skip:]
    nz = nz[skip:]
    zx = zx[skip:]
    zy = zy[skip:]
    zz = zz[skip:]
    sx = sx[skip:]
    sy = sy[skip:]
    sz = sz[skip:]

    # CREO LE STORIE MONTECARLO DEI MODULI DEI VETTORI
    M = np.sqrt(mx**2 + my**2 + mz**2)
    N = np.sqrt(nx**2 + ny**2 + nz**2)
    Z = np.sqrt(zx**2 + zy**2 + zz**2)
    S = np.sqrt(sx**2 + sy**2 + sz**2)


    # CREO GLI ARRAY PER SALVARE L'ERRORE
    # AGGIUNGERE ANCHE LE ALTRE OSSERVABILI
    err_ene_dens = np.zeros(0)
    err_spec_heat = np.zeros(0)

    err_mx = np.zeros(0)
    err_my = np.zeros(0)
    err_mz = np.zeros(0)
    err_nx = np.zeros(0)
    err_ny = np.zeros(0)
    err_nz = np.zeros(0)
    err_zx = np.zeros(0)
    err_zy = np.zeros(0)
    err_zz = np.zeros(0)
    err_sx = np.zeros(0)
    err_sy = np.zeros(0)
    err_sz = np.zeros(0)

    err_M = np.zeros(0)
    err_N = np.zeros(0)
    err_Z = np.zeros(0)
    err_S = np.zeros(0)

    err_binder_M = np.zeros(0)
    err_binder_N = np.zeros(0)
    err_binder_Z = np.zeros(0)
    err_binder_S = np.zeros(0)

    err_susc_M = np.zeros(0)
    err_susc_N = np.zeros(0)
    err_susc_Z = np.zeros(0)
    err_susc_S = np.zeros(0)

    # FACCIO DUE VOLTE LA PROCEDURA DI JACKKNIFE E POI MEDIO SU DUE TAGLIE DIVERSE
    taglie = np.array([len(ene_dens)//300, len(ene_dens)//100, len(ene_dens)//200, len(ene_dens)//500])

    for size in taglie:
        err_ene_dens = np.append(err_ene_dens, jack_error_1obs(media, ene_dens, size))
        err_spec_heat = np.append(err_spec_heat, jack_error_1obs(spec_heat, ene_dens, size))

        err_mx = np.append(err_mx, jack_error_1obs(media, mx, size))
        err_my = np.append(err_my, jack_error_1obs(media, my, size))
        err_mz = np.append(err_mz, jack_error_1obs(media, mz, size))
        err_nx = np.append(err_nx, jack_error_1obs(media, nx, size))
        err_ny = np.append(err_ny, jack_error_1obs(media, ny, size))
        err_nz = np.append(err_nz, jack_error_1obs(media, nz, size))
        err_zx = np.append(err_zx, jack_error_1obs(media, zx, size))
        err_zy = np.append(err_zy, jack_error_1obs(media, zy, size))
        err_zz = np.append(err_zz, jack_error_1obs(media, zz, size))
        err_sx = np.append(err_sx, jack_error_1obs(media, sx, size))
        err_sy = np.append(err_sy, jack_error_1obs(media, sy, size))
        err_sz = np.append(err_sz, jack_error_1obs(media, sz, size))

        err_M = np.append(err_M, jack_error_1obs(media, M, size))
        err_N = np.append(err_N, jack_error_1obs(media, N, size))
        err_Z = np.append(err_Z, jack_error_1obs(media, Z, size))
        err_S = np.append(err_S, jack_error_1obs(media, S, size))

        err_binder_M = np.append(err_binder_M, jack_error_1obs(binder, M, size))
        err_binder_N = np.append(err_binder_N, jack_error_1obs(binder, N, size))
        err_binder_Z = np.append(err_binder_Z, jack_error_1obs(binder, Z, size))
        err_binder_S = np.append(err_binder_S, jack_error_1obs(binder, S, size))

        err_susc_M = np.append(err_susc_M, jack_error_3obs(susc, mx, my, mz, size))
        err_susc_N = np.append(err_susc_N, jack_error_3obs(susc, nx, ny, nz, size))
        err_susc_Z = np.append(err_susc_Z, jack_error_3obs(susc, zx, zy, zz, size))
        err_susc_S = np.append(err_susc_S, jack_error_3obs(susc, sx, sy, sz, size))

    # STAMPO SU FILE I RISULTATI NEL FORMATO ENE_SP, ENE_G, ENE_DENS, SUSC, GPM, C, BINDER, CSI
    if(os.path.isdir('./data_w_errors') == True):
        f_name = "L_%d_alpha_%f" % (L, alpha)
        output_file = open('./data_w_errors/%s.dat' % (f_name) , 'a')
    else:
        dir_name = "./data_w_errors"
        os.makedirs(dir_name, exist_ok = False)
        f_name = "L_%d_alpha_%f" % (L, alpha)
        output_file = open('./data_w_errors/%s.dat' % (f_name) , 'a')

    output_file.write(str(L) + '\t')                        # COLONNA 1
    output_file.write(str(T) + '\t')                        # COLONNA 2
    output_file.write(str(alpha) + '\t')                    # COLONNA 3

    output_file.write(str(np.mean(ene_dens)) + '\t')        # COLONNA 4
    output_file.write(str(np.mean(err_ene_dens)) + '\t')    # COLONNA 5

    output_file.write(str(spec_heat(ene_dens)) + '\t')      # COLONNA 6
    output_file.write(str(np.mean(err_spec_heat)) + '\t')   # COLONNA 7

    output_file.write(str(np.mean(M)) + '\t')               # COLONNA 8
    output_file.write(str(np.mean(err_M)) + '\t')           # COLONNA 9
    output_file.write(str(np.mean(N)) + '\t')               # COLONNA 10
    output_file.write(str(np.mean(err_N)) + '\t')           # COLONNA 11
    output_file.write(str(np.mean(Z)) + '\t')               # COLONNA 12
    output_file.write(str(np.mean(err_Z)) + '\t')           # COLONNA 13
    output_file.write(str(np.mean(S)) + '\t')               # COLONNA 14
    output_file.write(str(np.mean(err_S)) + '\t')           # COLONNA 15

    output_file.write(str(binder(M)) + '\t')                # COLONNA 16
    output_file.write(str(np.mean(err_binder_M)) + '\t')    # COLONNA 17
    output_file.write(str(binder(N)) + '\t')                # COLONNA 18
    output_file.write(str(np.mean(err_binder_N)) + '\t')    # COLONNA 19
    output_file.write(str(binder(Z)) + '\t')                # COLONNA 20
    output_file.write(str(np.mean(err_binder_Z)) + '\t')    # COLONNA 21
    output_file.write(str(binder(S)) + '\t')                # COLONNA 22
    output_file.write(str(np.mean(err_binder_S)) + '\t')    # COLONNA 23

    output_file.write(str(susc(mx,my,mz)) + '\t')           # COLONNA 24
    output_file.write(str(np.mean(err_susc_M)) + '\t')      # COLONNA 25
    output_file.write(str(susc(nx,ny,nz)) + '\t')           # COLONNA 26
    output_file.write(str(np.mean(err_susc_N)) + '\t')      # COLONNA 27
    output_file.write(str(susc(zx,zy,zz)) + '\t')           # COLONNA 28
    output_file.write(str(np.mean(err_susc_Z)) + '\t')      # COLONNA 29
    output_file.write(str(susc(sx,sy,sz)) + '\t')           # COLONNA 30
    output_file.write(str(np.mean(err_susc_S)) + '\t')      # COLONNA 31

    output_file.write(str(np.mean(mx)) + '\t')              # COLONNA 32
    output_file.write(str(np.mean(err_mx)) + '\t')          # COLONNA 33
    output_file.write(str(np.mean(my)) + '\t')              # COLONNA 34
    output_file.write(str(np.mean(err_my)) + '\t')          # COLONNA 35
    output_file.write(str(np.mean(mz)) + '\t')              # COLONNA 36
    output_file.write(str(np.mean(err_mz)) + '\t')          # COLONNA 37
    output_file.write(str(np.mean(nx)) + '\t')              # COLONNA 38
    output_file.write(str(np.mean(err_nx)) + '\t')          # COLONNA 39
    output_file.write(str(np.mean(ny)) + '\t')              # COLONNA 40
    output_file.write(str(np.mean(err_ny)) + '\t')          # COLONNA 41
    output_file.write(str(np.mean(nz)) + '\t')              # COLONNA 42
    output_file.write(str(np.mean(err_nz)) + '\t')          # COLONNA 43
    output_file.write(str(np.mean(zx)) + '\t')              # COLONNA 44
    output_file.write(str(np.mean(err_zx)) + '\t')          # COLONNA 45
    output_file.write(str(np.mean(zy)) + '\t')              # COLONNA 46
    output_file.write(str(np.mean(err_zy)) + '\t')          # COLONNA 47
    output_file.write(str(np.mean(zz)) + '\t')              # COLONNA 48
    output_file.write(str(np.mean(err_zz)) + '\t')          # COLONNA 49
    output_file.write(str(np.mean(sx)) + '\t')              # COLONNA 50
    output_file.write(str(np.mean(err_sx)) + '\t')          # COLONNA 51
    output_file.write(str(np.mean(sy)) + '\t')              # COLONNA 52
    output_file.write(str(np.mean(err_sy)) + '\t')          # COLONNA 53
    output_file.write(str(np.mean(sz)) + '\t')              # COLONNA 54
    output_file.write(str(np.mean(err_sz)) + '\n')          # COLONNA 55

    output_file.close()

# CONTROLLO CON PLOT
# pl.xscale('log')
# pl.yscale('log')
# pl.scatter(taglie, err_ene_dens, label='ene_dens')
# pl.scatter(taglie, err_spec_heat, label='err_spec_heat')

# pl.scatter(taglie, err_mx, label ='mx')
# pl.scatter(taglie, err_my, label='my')
# pl.scatter(taglie, err_mz, label='mz')
# pl.scatter(taglie, err_nx, label='nx')
# pl.scatter(taglie, err_ny, label='ny')
# pl.scatter(taglie, err_nz, label='nz')
# pl.scatter(taglie, err_zx, label='zx')
# pl.scatter(taglie, err_zy, label='zy')
# pl.scatter(taglie, err_zz, label='zz')
# pl.scatter(taglie, err_sx, label='sx')
# pl.scatter(taglie, err_sy, label='sy')
# pl.scatter(taglie, err_sz, label='sz')

# pl.scatter(taglie, err_M, label='M')
# pl.scatter(taglie, err_N, label='N')
# pl.scatter(taglie, err_Z, label='Z')
# pl.scatter(taglie, err_S, label='S')

# pl.scatter(taglie, err_binder_M, label='err_binder_M')
# pl.scatter(taglie, err_binder_N, label='err_binder_N')
# pl.scatter(taglie, err_binder_Z, label='err_binder_Z')
# pl.scatter(taglie, err_binder_S, label='err_binder_S')

# pl.scatter(taglie, err_susc_M, label='err_susc_M')
# pl.scatter(taglie, err_susc_N, label='err_susc_N')
# pl.scatter(taglie, err_susc_Z, label='err_susc_Z')
# pl.scatter(taglie, err_susc_S, label='err_susc_S')
# pl.scatter(taglie, err_G_pm)
# pl.scatter(taglie, err_spec_heat)
# pl.scatter(taglie, err_binder)
# pl.scatter(taglie, err_corr_len)
# pl.scatter(taglie, err_ene_g)
# print('susc_N = ', susc(nx,ny,nz), '+-', np.mean(err_susc_N))
# pl.legend()
# pl.show()
