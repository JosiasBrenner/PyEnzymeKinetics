from numpy import max, where, abs, argmin, mean


def get_v(substrate_conc, time):
    v_all = 0.0*substrate_conc[:]  # initialize velocity vector
    if len(substrate_conc.shape) > 1:
        for i in range(substrate_conc.shape[0]):

            prev_value = substrate_conc[i, 0]
            prev_time = 0.0

            for j in range(substrate_conc.shape[1]):

                if time[j] == 0:
                    delta = prev_value - substrate_conc[i, j]
                else:
                    delta = abs(
                        (prev_value - substrate_conc[i, j])/(time[j]-prev_time))

                v_all[i, j] = delta
                prev_value = substrate_conc[i, j]
                prev_time = time[j]

        v = max(v_all, axis=0)

    else:

        prev_value = substrate_conc[0]
        prev_time = 0.0

        for j in range(substrate_conc.shape[0]):

            if time[j] == 0:
                delta = prev_value - substrate_conc[j]
            else:
                delta = abs(
                    (prev_value - substrate_conc[j])/(time[j]-prev_time))

            v_all[j] = delta
            prev_value = substrate_conc[j]
            prev_time = time[j]

        v = v_all
        print("done")

    return v


def get_initial_vmax(substrate_conc, time):
    v = get_v(substrate_conc, time)
    return max(v)


def get_initial_Km(substrate_conc, time):

    v = get_v(substrate_conc, time)
    idx_max = where(v == max(v))[0][0]
    idx_Km = (abs(v[idx_max:]-max(v)/2)).argmin()

    if len(substrate_conc.shape) > 1:
        km = mean(substrate_conc, axis=0)[idx_max+idx_Km]
    else:
        km = substrate_conc[idx_max+idx_Km]

    return km
