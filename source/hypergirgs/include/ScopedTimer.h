/**
 * @file
 * @brief ScopedTimer
 *
 * Simple timer to measure wall-time. By default its
 * measures the time span between creation and deletion.
 *
 * @author Manuel Penschuck
 * @copyright
 * Copyright (C) 2017 Manuel Penschuck
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * @copyright
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * @copyright
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once
#ifndef SCOPED_TIMER_HPP
#define SCOPED_TIMER_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <chrono>

class ScopedTimer {
    using Clock = std::chrono::high_resolution_clock;

    std::string _prefix;
    Clock::time_point _begin;

    uint64_t _scale;
    double _offset;

    double* _output;

public:
    ScopedTimer()
          : _begin(Clock::now()), _scale(1), _offset(0), _output(nullptr)
    {}

    ScopedTimer(const std::string& prefix, uint64_t scale=0, double offset=0.0)
          : _prefix(prefix), _begin(Clock::now()), _scale(scale), _offset(offset), _output(nullptr)
    {}

    ScopedTimer(double& output) : _begin(Clock::now()), _scale(1), _offset(0), _output(&output) {}

    ~ScopedTimer() {
        if (!_prefix.empty())
            report();

        if (_output)
            *_output = elapsed();
    }

    void start() {
        _begin = Clock::now();
    }

    double elapsed() const {
        const auto t2 = Clock::now();
        std::chrono::duration<double> time_span =
                std::chrono::duration_cast<std::chrono::duration<double>>(t2 - _begin);

        return (time_span.count()*1e3) - _offset;
    }

    double report() const {
        return report(_prefix);
    }

    double report(const std::string & prefix) const {
        const double timeUs = elapsed();
        
        if (!_scale) {
            std::cout << prefix << " Time elapsed: " << timeUs << "ms" << std::endl;
        } else {
            std::cout << prefix << " Time elapsed: " << timeUs << "ms / " << _scale << " = " <<  (1e3*timeUs / _scale) << "us" << std::endl;
        }

        return timeUs;
    }
};

#endif
