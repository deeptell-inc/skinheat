"""Consistency check (15% tolerance: mesoscopic leads have finite band W=4 and level-spacing broadening; the GF uses wide-band leads)."""
import numpy as np
from skinheat import build_chain_with_leads, ness_covariance, lead_currents, chain_heff, delta_heat


def lyap_dQ(L, kappa, phi, TL, TR, T0, K=160, W=4.0):
    def QR(TL_, TR_):
        m = build_chain_with_leads(L, kappa=kappa, phi=phi, K=K, W=W, TL=TL_, TR=TR_)
        c = lead_currents(m, ness_covariance(m))
        return -c["R"]["E"], -c["L"]["E"]
    QRf, _ = QR(TL, TR); QR0, _ = QR(T0, T0)
    _, QLr = QR(TR, TL); _, QL0 = QR(T0, T0)
    return QRf - QR0, QLr - QL0


def test_consistency():
    L, kappa, phi, Th, Tc = 12, 0.3, 0.7, 0.8, 0.5
    dQf, dQr = lyap_dQ(L, kappa, phi, Th, Tc, Tc)
    gf, gr, *_ = delta_heat(chain_heff(L, kappa=kappa, phi=phi), Th, Tc)
    assert abs(dQf - gf) < 0.15 * max(abs(gf), 1e-6), (dQf, gf)
    assert abs(dQr - gr) < 0.15 * max(abs(gr), 1e-6), (dQr, gr)


if __name__ == "__main__":
    test_consistency(); print("OK")
