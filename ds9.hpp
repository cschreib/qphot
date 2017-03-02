bool read_ds9_region_circles(std::string file_name, vec2d& regs, bool& physical,
    std::string color) {
    if (!file::exists(file_name)) {
        write_error("could not open region file '", file_name, "'");
        return false;
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
            write_error(file_name, ":", l, ": ",
                "ill formed 'circle' line, expecting 3 arguments, got ", args.size());
            return false;
        }

        double ra, dec, rad;
        args = trim(args);
        if (args[0].find_first_of(':') != args[0].npos) {
            if (!sex2deg(args[0], args[1], ra, dec)) {
                write_error(file_name, ":", l, ": ",
                    "could not convert sexagesimal coordinates to degrees");
                return false;
            }

            if (!end_with(args[2], "\"")) {
                write_error(file_name, ":", l, ": expected radius in arcsec");
                return false;
            }
        } else {
            if (!from_string(args[0], ra) || !from_string(args[1], dec)) {
                write_error(file_name, ":", l, ": ",
                    "could not read coordinates to ", (physical ? "(x,y)" : "degrees"));
                return false;
            }
        }

        if (physical) {
            if (!from_string(args[2], rad)) {
                write_error(file_name, ":", l, ": could not read radius in pixels");
                return false;
            }
        } else {
            args[2] = erase_end(args[2], "\"");
            if (!from_string(args[2], rad)) {
                write_error(file_name, ":", l, ": could not read radius in arcsec");
                return false;
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

    return true;
}

bool read_ds9_region_circles_physical(std::string file_name,
    const astro::wcs& w, vec2d& regs, std::string color = "") {

    bool physical = false;
    if (!read_ds9_region_circles(file_name, regs, physical, color)) return false;

    if (physical) return true;

    if (!w.is_valid()) {
        return false;
    }

    for (uint_t i : range(regs.dims[0])) {
        double x, y;
        astro::ad2xy(w, regs(i,0), regs(i,1), x, y);
        regs(i,0) = x-1;
        regs(i,1) = y-1;
    }

    double aspix = 1.0;
    if (!astro::get_pixel_size(w, aspix)) {
        return false;
    }

    regs(_,2) /= aspix;

    return true;
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
