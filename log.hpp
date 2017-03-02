std::ofstream errfile;

template<typename T>
void write_impl(const T& msg) {
    errfile << msg;
}

template<typename T>
void swallow_arguments(std::initializer_list<T>) {}

template<typename ... Args>
void write_errfile(const Args& ... args) {
    if (!errfile.is_open()) return;

    // Print time stamp
    std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm tm = *std::localtime(&t);
    errfile << "["
            << align_right(strn(tm.tm_year+1900),4,'0') << "/"
            << align_right(strn(tm.tm_mon+1),2,'0') << "/"
            << align_right(strn(tm.tm_mday),2,'0') << " "
            << align_right(strn(tm.tm_hour),2,'0') << ":"
            << align_right(strn(tm.tm_min),2,'0') << ":"
            << align_right(strn(tm.tm_sec),2,'0') << "]: ";

    // Print arguments
    swallow_arguments({(write_impl(args), 0)...});

    // End line
    errfile << std::endl;
}

template<typename ... Args>
void write_error(const Args& ... args) {
    error(args...);
    write_errfile("error: ", args...);
}

template<typename ... Args>
void write_warning(const Args& ... args) {
    warning(args...);
    write_errfile("warning: ", args...);
}
