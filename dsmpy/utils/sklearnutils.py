from dsmpy.seismicmodel import SeismicModel
from dsmpy.modelparameters import ModelParameters, ParameterType
from dsmpy.windowmaker import WindowMaker
from dsmpy.component import Component
from dsmpy.dataset import Dataset
from dsmpy.dsm import compute_models_parallel
from dsmpy import root_sac
from pytomo.work.ca.params import get_dataset, get_model_syntest1_prem
from pytomo.preproc.iterstack import find_best_shift
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from mpi4py import MPI
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import logging

# logging.basicConfig(
#     level=logging.INFO, filename='sklearnutils.log', filemode='w')

def get_XY(
        model, dataset, windows, tlen, nspc,
        freq, freq2, filter_type='bandpass',
        sampling_hz=5, var=2.5, ratio=2.5, corr=0.,
        phase_ref=None, buffer=10.,
        mode=0) -> (np.ndarray, np.ndarray):
    """Compute the feature matrix X and target vector y to be used
    as input to scikit-learn linear models.

    X and y are linked by the equation Xm = y, where m is the model
    parameter vector. The order for m is given by the order in
    SeismicModel.gradient_models(), and is:
    [[radial_nodes_for_type_1] + [radial_nodes_for_type_2] + ...].

    This method should be able to scale to large dataset, since
    the computations are done in the frequency domain (typically approx.
    256 to 512 np.Complex64 per synthetic), the transformation to time
    domain is done event by event (the data is freed after), and the
    gradient matrix X contains windowed time series with typically
    a few hundreth to thousands of floats. Furthermore, only the
    frequency-domain synthetics are replicated on all cores. All the
    time domain operations, as well as X and y are defined on thread 0
    only. For instance, 10,000 records sampled at 5 Hz for 50 s windows
    for one seismic component with 100 model parameters should not take
    more than approx. 1e4 * 5 * 50 * 101 * 6.4e-8 = 16.2 Gb.

    Args:
        model (SeismicModel): model at which the gradient is evaluated.
            Must be a mesh and have model._model_params not None.
        dataset (Dataset): dataset
        windows (list of Window): time windows
        tlen (float): length of time series for synthetics
        nspc (int): number of points in frequency domain for synthetics
        sampling_hz (int): sampling frequency of synthetics in time
            domain. Better to divide 20.
        var (float): variance cutoff. Records with variance > var
            will be excluded (default is 2.5).
        ratio (float): amplitude ratio cutoff. Records with
            1/(obs/syn) < ratio or obs/syn > ratio will be excluded
            (default is 2.5).
        corr (float): correlation coefficient cutoff. Records with
            correlation < corr will be excluded (default is 0).
        phase_ref (str): reference phase for static correction
            (default is None).
        buffer (float): time buffer in seconds for static
            correction (default is 10).
        mode (int): commputation mode. 0: P-SV + SH, 1: P-SV, 2: SH

    Returns:
        np.ndarray: X, the waveform gradient with respec to model.
            The shape is (n_time_points, n_model_parameters).
        np.ndarray: y, the waveform residual vector.
            The shape is (n_time_points,).
    """
    if model.get_model_params() is None:
        return None, None
    grad_models, dxs = model.gradient_models()
    outputs = compute_models_parallel(
        dataset, grad_models + [model], tlen=tlen,
        nspc=nspc, sampling_hz=dataset.sampling_hz, mode=mode)

    sample_skip = int(dataset.sampling_hz // sampling_hz)

    if MPI.COMM_WORLD.Get_rank() == 0:
        grad_outputs = outputs[:len(grad_models)]
        ref_outputs = outputs[len(grad_models)]

        n_params = len(grad_models) // 2
        n_ev = len(dataset.events)

        buffer_ = int(buffer * dataset.sampling_hz)
        shift_dict = dict()
        checked = set()

        logging.info('Building y')
        y = []
        for iev in range(n_ev):
            event = dataset.events[iev]
            ref_output = ref_outputs[iev]
            ref_output.filter(freq, freq2, filter_type)
            start, end = dataset.get_bounds_from_event_index(iev)
            for ista in range(start, end):
                station = dataset.stations[ista]
                jsta = np.argwhere(ref_output.stations == station)[0][0]
                windows_filt = [
                    window for window in windows
                    if (window.station == station
                        and window.event == event)]

                # shift window to maximize correlation
                if phase_ref is not None:
                    for iwin, window in enumerate(windows_filt):
                        if window.phase_name == phase_ref:
                            window_arr = window.to_array()
                            icomp = window.component.value
                            i_start = int(window_arr[0] * dataset.sampling_hz)
                            i_end = int(window_arr[1] * dataset.sampling_hz)
                            data_cut = dataset.data[
                                       iwin, icomp, ista, :i_end - i_start]
                            u_cut_buffer = ref_output.us[
                               icomp, jsta, i_start + buffer_:i_end - buffer_]
                            shift, _ = find_best_shift(
                                data_cut, u_cut_buffer, skip_freq=2)
                            key = _get_shift_key(window, phase_ref)
                            shift_dict[key] = shift

                            # window_arr = window.to_array()
                            # i_start_syn = int(
                            #     window_arr[0] * dataset.sampling_hz) + buffer_
                            # i_end_syn = int(
                            #     window_arr[1] * dataset.sampling_hz) - buffer_
                            # i_start_dat = shift_dict[key]
                            # i_end_dat = i_start_dat + i_end_syn - i_start_syn
                            # print(i_start_dat, i_end_dat)
                            # u_cut = ref_output.us[
                            #     icomp, jsta, i_start_syn:i_end_syn]
                            # data_cut = dataset.data[
                            #     iwin, icomp, ista, i_start_dat:i_end_dat]
                            # plt.plot(u_cut, label='ref')
                            # plt.plot(data_cut)
                            # plt.legend()
                            # plt.show()

                for iwin, window in enumerate(windows_filt):
                    if window.phase_name == phase_ref:
                        continue
                    key = _get_shift_key(window, phase_ref)
                    if phase_ref is not None and key not in shift_dict:
                        continue

                    window_arr = window.to_array()
                    i_start_syn = int(
                        window_arr[0] * dataset.sampling_hz) + buffer_
                    i_end_syn = int(
                        window_arr[1] * dataset.sampling_hz) - buffer_
                    i_start_dat = shift_dict[key]
                    i_end_dat = i_start_dat + i_end_syn - i_start_syn
                    ref_u = ref_output.us[
                        icomp, jsta, i_start_syn:i_end_syn]
                    data_cut = dataset.data[
                        iwin, icomp, ista, i_start_dat:i_end_dat]

                    # plt.plot(ref_u, label='ref')
                    # plt.plot(data_cut)
                    # plt.legend()
                    # plt.show()

                    if check_data(data_cut, ref_u, var, corr, ratio):
                        key = get_check_key(window)
                        checked.add(key)
                        residual = data_cut - ref_u
                        residual /= np.abs(data_cut).max()
                        y.append(residual[::sample_skip])

            ref_output.free()

        logging.info('Building X')
        _, itypes, iradii = model.get_model_params().get_free_all_indices()
        types = [model.get_model_params().get_types()[i]
                 for i in itypes]
        radii = [model.get_model_params().get_grd_params()[i]
                 for i in iradii]
        x = []
        for imod in range(n_params):
            logging.info(f'model {imod}')
            x_imod = []
            for iev in range(n_ev):
                event = dataset.events[iev]
                output_plus = grad_outputs[imod][iev]
                output_minus = grad_outputs[imod + n_params][iev]
                output_plus.filter(freq, freq2, filter_type)
                output_minus.filter(freq, freq2, filter_type)
                start, end = dataset.get_bounds_from_event_index(iev)
                for ista in range(start, end):
                    station = dataset.stations[ista]
                    jsta = np.argwhere(output_plus.stations == station)[0][0]
                    windows_filt = [
                        window for window in windows
                        if (window.station == station
                            and window.event == event)]

                    for iwin, window in enumerate(windows_filt):
                        if window.phase_name == phase_ref:
                            continue
                        key = _get_shift_key(window, phase_ref)
                        if phase_ref is not None and key not in shift_dict:
                            continue

                        icomp = window.component.value
                        window_arr = window.to_array()
                        i_start_dat = shift_dict[key]
                        i_end_dat = i_start_dat + i_end_syn - i_start_syn
                        i_start_syn = int(
                            window_arr[0] * dataset.sampling_hz) + buffer_
                        i_end_syn = int(
                            window_arr[1] * dataset.sampling_hz) - buffer_
                        data_cut = dataset.data[
                                   iwin, icomp, ista, i_start_dat:i_end_dat]

                        check_key = get_check_key(window)
                        if check_key in checked:
                            grad_u = (
                                output_plus.us[
                                    icomp, jsta, i_start_syn:i_end_syn]
                                - output_minus.us[
                                    icomp, jsta, i_start_syn:i_end_syn]
                            ) / (dxs[imod] - dxs[imod + n_params])
                            grad_u /= np.abs(data_cut).max()

                            # if QMU or QKAPPA, compute gradient
                            # with respect to 1/Q
                            if types[imod] == ParameterType.QMU:
                                r = radii[imod]
                                Qmu = model.get_value_at(r, ParameterType.QMU)
                                grad_u *= -Qmu * Qmu
                            if types[imod] == ParameterType.QKAPPA:
                                r = radii[imod]
                                Qk = model.get_value_at(r, ParameterType.QKAPPA)
                                grad_u *= -Qk * Qk

                            x_imod.append(grad_u[::sample_skip])

                output_plus.free()
                output_minus.free()
            x.append(np.hstack(x_imod))

        x = np.array(x).T
        y = np.hstack(y).reshape((-1, 1))

        return x, y, checked
    else:
        return None, None, None


def _get_shift_key(window, phase_ref):
    return str(window.station) + window.event.event_id + phase_ref

def get_check_key(window):
    return str(window.station) + window.event.event_id + window.phase_name + str(window.component)

def misfits(data, syn):
    """Returns variance, corr, ratio.
    """
    residual = data - syn
    mean = np.mean(residual)
    variance = np.dot(residual - mean, residual - mean)
    if np.abs(data).max() > 0:
        variance /= np.dot(data - data.mean(),
                           data - data.mean())
    else:
        variance = np.nan
    corr = np.corrcoef(data, syn)[0, 1]
    if syn.max() > 0:
        ratio = (data.max() - data.min()) / (
                syn.max() - syn.min())
    else:
        ratio = np.nan
    return variance, corr, ratio


def check_data(data, syn, var, corr, ratio):
    var_, corr_, ratio_ = misfits(data, syn)
    return (
            np.isfinite(var_)
            and np.isfinite(corr_)
            and np.isfinite(ratio_)
            and var_ <= var
            and corr_ >= corr
            and 1 / ratio <= ratio_ <= ratio
    )

if __name__ == '__main__':
    types = [ParameterType.VSH]
    radii = [3480. + 20 * i for i in range(41)]
    model_params = ModelParameters(
        types=types,
        radii=radii,
        mesh_type='boxcar')
    model = SeismicModel.prem().boxcar_mesh(model_params)

    n_param = model_params.get_n_params()
    its, itypes, igrds = model_params.get_free_all_indices()
    order_types = [model_params.get_types()[itypes[i]]
                   for i in its]
    order_radii = [model_params.get_grd_params()[igrds[i]]
                   for i in its]

    # sac_files = glob.glob(os.path.join(root_sac, '*[RZT]'))
    # dataset = Dataset.dataset_from_sac(
    #     sac_files, headonly=False)

    t_before = 20.
    t_after = 40.
    sampling_hz = 20
    window_npts = int((t_before + t_after) * 20)
    tlen = 1638.4
    nspc = 256
    mode = 2
    freq = 0.01
    freq2 = 0.08
    filter_type = 'bandpass'

    dataset, output = get_dataset(
        get_model_syntest1_prem(), tlen=tlen, nspc=nspc,
        sampling_hz=sampling_hz,
        mode=mode, add_noise=False, noise_normalized_std=1.)

    windows = WindowMaker.windows_from_dataset(
        dataset, 'prem', ['ScS'], [Component.T],
        t_before=t_before, t_after=t_after)

    dataset.filter(freq, freq2, filter_type)
    dataset.apply_windows(windows, 1, window_npts)

    X, y = get_XY(
        model, dataset, windows, tlen=tlen, nspc=nspc,
        freq=freq, freq2=freq2, filter_type=filter_type,
        sampling_hz=sampling_hz, mode=mode)

    if MPI.COMM_WORLD.Get_rank() == 0:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        n_pc = 6
        pca = PCA(n_components=n_pc)
        pca.fit(X_scaled)

        fig10, ax10 = plt.subplots()
        ax10.bar(list(range(n_pc)), pca.explained_variance_ratio_)

        fig11, ax11 = plt.subplots(6, 2)
        for i, ax in enumerate(ax11.ravel()):
            ax.bar(list(range(X.shape[1])), pca.components_[i])
        plt.show()

        X_scaled_pca = pca.transform(X_scaled)

        reg = linear_model.Ridge(alpha=0.)

        reg.fit(X_scaled_pca, y)
        y_pred = reg.predict(X_scaled_pca)
        mse = mean_squared_error(y, y_pred)
        rmse = np.sqrt(mse)
        print(rmse)

        coefs_scaled_pca = reg.coef_.reshape(1, -1)
        coefs_scaled = pca.inverse_transform(coefs_scaled_pca)
        coefs = scaler.transform(coefs_scaled)

        best_params = np.array(coefs).reshape(len(types), -1)
        value_dict = {p_type: best_params[i]
                      for i, p_type in enumerate(types)}
        values = model_params.get_values_matrix(value_dict)
        best_model = model.multiply(values)

        fig0, ax0 = plt.subplots(1)
        model.plot(types=types, ax=ax0, label='prem')
        get_model_syntest1_prem().plot(types=types, ax=ax0, label='target')
        best_model.plot(types=types, ax=ax0, label='inverted')
        ax0.set_ylim([3480, 4000])
        ax0.set_xlim([6.5, 8])
        ax0.legend()

        fig, ax = plt.subplots(1)
        ax.imshow(X, aspect='auto')

        fig2, ax2 = plt.subplots(1)
        ax2.imshow(np.dot(X.T, X))

        plt.show()

