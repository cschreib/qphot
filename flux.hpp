void get_background(const vec2d& data, double& bg, double& rms) {
    vec1u idg;
    idg = where(is_finite(data));
    bg = median(data[idg]);
    rms = 1.48*mad(data[idg]);

    idg = where(is_finite(data) && data - bg < 10*rms);
    bg = median(data[idg]);
    rms = 1.48*mad(data[idg]);
}

void get_background_apertures(const vec2b& mask, double radius, vec1u& y, vec1u& x) {
    y.clear();
    x.clear();

    // First segment the mask
    uint_t nseg = 0;
    vec2u seg = segment(mask, nseg);

    // Add segments at the edges
    seg(0,_) = seg(seg.dims[0]-1,_) = seg(_,0) = seg(_,seg.dims[1]-1) = nseg+1;

    // Build distance from segments
    vec2d dist;
    vec2u dist_id;
    segment_distance(seg, dist, dist_id);

    // Now add circles on the largest distance point repeatedly until there is no
    // space left
    double radius2 = sqr(radius);
    uint_t mid = max_id(dist);
    while (dist[mid] >= radius2) {
        vec1u muid = mult_ids(dist, mid);
        y.push_back(muid[0]);
        x.push_back(muid[1]);
        dist[mid] = 0;
        dist_id[mid] = mid;

        std::vector<uint_t> oy = {muid[0]}, ox = {muid[1]};
        while (!ox.empty()) {
            auto toy = oy, tox = ox;
            oy.clear(); ox.clear();

            while (!tox.empty()) {
                uint_t ty = toy.back(); toy.pop_back();
                uint_t tx = tox.back(); tox.pop_back();

                auto check_add = [&](uint_t tty, uint_t ttx) {
                    if (dist.safe(tty,ttx) == 0) return;

                    if (dist_id.safe(ty,tx) == mid) {
                        double nd = sqr(double(tty) - y.back())
                                  + sqr(double(ttx) - x.back());

                        if (nd < radius2) {
                            // Still within the circle, move on
                            dist.safe(tty,ttx) = 0;
                            dist_id.safe(tty,ttx) = mid;
                            oy.push_back(tty);
                            ox.push_back(ttx);
                        } else {
                            // Setup new origin and continue below
                            dist_id.safe(ty,tx) = flat_id(dist, ty, tx);
                        }
                    }

                    if (dist_id.safe(ty,tx) != mid) {
                        double x0 = dist_id.safe(ty,tx) % dist_id.dims[1];
                        double y0 = dist_id.safe(ty,tx) / dist_id.dims[1];
                        double nd = sqr(double(tty) - y0) + sqr(double(ttx) - x0);
                        if (dist.safe(tty,ttx) > nd) {
                            dist.safe(tty,ttx) = nd;
                            dist_id.safe(tty,ttx) = dist_id.safe(ty,tx);
                            oy.push_back(tty);
                            ox.push_back(ttx);
                        }
                    }
                };

                if (ty != 0)              check_add(ty-1,tx);
                if (ty != dist.dims[0]-1) check_add(ty+1,tx);
                if (tx != 0)              check_add(ty,tx-1);
                if (tx != dist.dims[1]-1) check_add(ty,tx+1);
            }
        }

        mid = max_id(dist);
    }
}

// Measure the "seeing" (FWHM) from a PSF image
double get_seeing(vec2d img, double y0, double x0) {
    auto get_seeing_raw = [](const vec2d& timg, double ty0, double tx0) {
        double tot = 0.0;
        for (uint_t y : range(timg.dims[0]))
        for (uint_t x : range(timg.dims[1])) {
            tot += timg.safe(y,x)*(sqr(y - ty0) + sqr(x - tx0));
        }

        return sqrt(tot);
    };

    // Compute first estimate
    img /= total(img); // image needs to be normalized first
    double rr = get_seeing_raw(img, y0, x0);

    // Simulate various profiles
    vec1d td = rgen_step(0.1, max(img.dims[0], img.dims[1]), 0.3);
    vec1d md(td.dims);

    vec2d timg(img.dims);
    for (uint_t i : range(td)) {
        for (uint_t y : range(timg.dims[0]))
        for (uint_t x : range(timg.dims[1])) {
            timg.safe(y,x) = integrate_gauss_2d(y-0.5, y+0.5, x-0.5, x+0.5, y0, x0, td[i]);
        }

        md[i] = get_seeing_raw(timg, y0, x0);
    }

    // Correct measured value using simulation
    // The additional factor of 2*sqrt(log(4.0)) is the conversion
    // from a Gaussian sigma into a FWHM (=2.355)
    return interpolate(td, md, rr)*2*sqrt(log(4.0));
}

vec2d make_point_source_integer(const vec2d& psf, const std::array<uint_t,2>& dims, int_t y, int_t x) {
    vec2d model(dims);
    int_t y0 = max(0, y - int_t(psf.dims[0]/2));
    int_t y1 = min(int_t(dims[0]), y + int_t(psf.dims[0]/2+1));
    int_t x0 = max(0, x - int_t(psf.dims[1]/2));
    int_t x1 = min(int_t(dims[0]), x + int_t(psf.dims[1]/2+1));

    for (int_t iy = y0; iy < y1; ++iy)
    for (int_t ix = x0; ix < x1; ++ix) {
        int_t ipy = iy - y + int_t(psf.dims[0]/2);
        int_t ipx = ix - x + int_t(psf.dims[1]/2);
        model.safe(uint_t(iy),uint_t(ix)) = psf.safe(uint_t(ipy),uint_t(ipx));
    }

    return model;
}

vec2d make_point_source(const vec2d& psf, const std::array<uint_t,2>& dims, double y, double x) {
    int_t iy = round(y), ix = round(x);
    vec2d npsf = translate(psf, y - iy, x - ix);
    return make_point_source_integer(npsf, dims, iy, ix);
}
