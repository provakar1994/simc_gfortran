void boost_ep_to_breit() {
    double M_p = 0.938; // Proton mass

    // Beam electron: 10 GeV along z
    TLorentzVector k_in(0, 0, 10.0, 10.0);

// Scattered electron: 8 GeV at 20 deg
double E_out = 8.0;
double theta = TMath::DegToRad() * 20.0;
    double px = E_out * sin(theta);
    double pz = E_out * cos(theta);
    TLorentzVector k_out(px, 0, pz, E_out);

    // Proton at rest
    TLorentzVector p(0, 0, 0, M_p);

    // q = k_in - k_out
    TLorentzVector q = k_in - k_out;
    double Q2 = -q.Mag2();
    double x_bj = Q2 / (2 * p.Dot(q));

    std::cout << "Q^2 = " << Q2 << " GeV^2" << std::endl;
    std::cout << "Bjorken x = " << x_bj << std::endl;

    // Define Breit frame boost direction (along q)
    TVector3 qvec = q.Vect();
    double q0 = q.E();

    // We want to boost to frame where q0 = 0
    // In such frame, the proton will move with -qvec/2
    // Compute boost vector to Breit frame:
    TVector3 beta_breit = qvec * (1.0 / q0); // β = q⃗ / q⁰

    // Check magnitude
    if (beta_breit.Mag() >= 1.0) {
        std::cerr << "Invalid boost: beta >= 1!" << std::endl;
        return;
    }

    // Boost all four-vectors
    TLorentzVector k_in_b = k_in;
    TLorentzVector k_out_b = k_out;
    TLorentzVector p_b = p;
    TLorentzVector q_b = q;

    // Boost to Breit frame (reverse direction)
    k_in_b.Boost(-beta_breit);
    k_out_b.Boost(-beta_breit);
    p_b.Boost(-beta_breit);
    q_b.Boost(-beta_breit);

    std::cout << "\n--- In Breit Frame ---" << std::endl;
    std::cout << "q0 = " << q_b.E() << " (should be ~0)" << std::endl;
    std::cout << "qz = " << q_b.Pz() << std::endl;
    std::cout << "Incoming electron: "; k_in_b.Print();
    std::cout << "Scattered electron: "; k_out_b.Print();
    std::cout << "Proton: "; p_b.Print();
}
