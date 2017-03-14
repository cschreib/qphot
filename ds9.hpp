void read_ds9_region_circles(std::string file_name, vec2d& regs, bool& physical,
    std::string color) {

    if (!file::exists(file_name)) {
        throw std::runtime_error("could not open region file '"+file_name+"'");
    }

    std::ifstream file(file_name);

    std::string global_color = "green";
    physical = false;
    std::string line;
    uint_t l = 0;
    while (std::getline(file, line)) {
        ++l;
        if (line.empty() || trim(line).empty() || trim(line)[0] == '#') continue;

        if (start_with(line, "global")) {
            std::string key = "color=";
            auto pos = line.find(key);
            if (pos != line.npos) {
                pos += key.size();
                global_color = trim(line.substr(pos, line.find_first_of(" \t", pos)-pos));
            }
            continue;
        }

        auto spos = line.find_first_of('(');
        if (spos == line.npos) {
            if (trim(line) == "fk5") physical = false;
            if (trim(line) == "physical") physical = true;
            continue;
        }

        std::string type = trim(line.substr(0, spos));
        if (type != "circle") continue;

        auto epos = line.find_first_of(')', spos+1);
        std::string targs = line.substr(spos+1, epos-(spos+1));
        vec1s args = split(targs, ",");
        if (args.size() != 3) {
            throw std::runtime_error(file_name+":"+strn(l)+": "
                "ill formed 'circle' line, expecting 3 arguments, got "+strn(args.size()));
        }

        double ra, dec, rad;
        args = trim(args);
        if (args[0].find_first_of(':') != args[0].npos) {
            if (!sex2deg(args[0], args[1], ra, dec)) {
                throw std::runtime_error(file_name+":"+strn(l)+": "
                    "could not convert sexagesimal coordinates to degrees");
            }

            if (!end_with(args[2], "\"")) {
                throw std::runtime_error(file_name+":"+strn(l)+": expected radius in arcsec");
            }
        } else {
            if (!from_string(args[0], ra) || !from_string(args[1], dec)) {
                throw std::runtime_error(file_name+":"+strn(l)+": "
                    "could not read coordinates to "+std::string(physical ? "(x,y)" : "degrees"));
            }
        }

        if (physical) {
            if (!from_string(args[2], rad)) {
                throw std::runtime_error(file_name+":"+strn(l)+": could not read radius in pixels");
            }
        } else {
            args[2] = erase_end(args[2], "\"");
            if (!from_string(args[2], rad)) {
                throw std::runtime_error(file_name+":"+strn(l)+": could not read radius in arcsec");
            }
        }

        if (!color.empty()) {
            std::string rcol = global_color;
            spos = line.find_first_of('#', epos);
            if (spos != line.npos) {
                std::string key = "color=";
                auto pos = line.find(key, spos+1);
                if (pos != line.npos) {
                    pos += key.size();
                    rcol = trim(line.substr(pos, line.find_first_of(" \t", pos)-pos));
                }
            }

            if (rcol != color) continue;
        }

        append<0>(regs, vec2d{{ra, dec, rad}});
    }
}

void read_ds9_region_circles_physical(std::string file_name,
    const astro::wcs& w, vec2d& regs, std::string color = "") {

    bool physical = false;
    read_ds9_region_circles(file_name, regs, physical, color);

    if (physical) return;

    double aspix = 1.0;
    if (!w.is_valid() || !astro::get_pixel_size(w, aspix)) {
        throw std::runtime_error("invalid WCS, cannot convert regions");
    }

    for (uint_t i : range(regs.dims[0])) {
        double x, y;
        astro::ad2xy(w, regs(i,0), regs(i,1), x, y);
        regs(i,0) = x-1;
        regs(i,1) = y-1;
    }

    regs(_,2) /= aspix;
}

void write_ds9_region(const std::string& filename, const astro::wcs& w,
    const vec1u& y, const vec1u& x, double radius) {

    vec1d ra, dec;
    astro::xy2ad(w, x+1.0, y+1.0, ra, dec);
    double aspix;
    astro::get_pixel_size(w, aspix);
    radius *= aspix;

    std::ofstream regfile(filename);
    regfile << "# Region file format: DS9 version 4.1\n";
    regfile << "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal "
        "roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 "
        "source=1\n";
    regfile << "fk5\n";

    regfile.precision(12);
    for (uint_t i : range(x)) {
        regfile << "circle(" << ra[i] << ", " << dec[i] << ", " << radius << "\")\n";
    }
}
