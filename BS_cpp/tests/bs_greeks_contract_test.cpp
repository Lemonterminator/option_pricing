#include "bs_batch.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace black_scholes;

namespace {
struct Row {
    std::string name;
    std::string kind;
    std::array<double, 6> inputs;
    ScalarGreeks expected;
};

double number(std::istream& stream) {
    std::string field;
    if (!std::getline(stream, field, ',')) throw std::runtime_error("Missing CSV number");
    char* end = nullptr;
    const double value = std::strtod(field.c_str(), &end);
    if (end == field.c_str() || (*end && *end != '\r')) throw std::runtime_error("Invalid CSV number");
    return value;
}

std::vector<Row> read_cases(const char* path) {
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("Cannot read fixture");
    std::string line;
    std::getline(stream, line);  // Header
    std::vector<Row> rows;
    while (std::getline(stream, line)) {
        std::istringstream fields(line);
        Row row;
        std::getline(fields, row.name, ',');
        std::getline(fields, row.kind, ',');
        if (row.kind != "call" && row.kind != "put") throw std::runtime_error("Invalid kind");
        for (double& input : row.inputs) input = number(fields);
        row.expected = {number(fields), number(fields), number(fields)};
        rows.push_back(row);
    }
    if (rows.empty()) throw std::runtime_error("No fixture rows");
    return rows;
}

bool close(double actual, double expected) {
    if (actual == expected) return true;
    if (!(std::isfinite(actual) && std::isfinite(expected))) return false;
    return std::abs(actual - expected) <= 2e-10 * std::abs(expected)
        + 4 * std::numeric_limits<double>::denorm_min();
}

void check(const Row& row, const char* method, ScalarGreeks actual, bool dump) {
    if (dump) {
        std::cout << row.name << ',' << row.kind << ',' << method << ','
                  << actual.vega << ',' << actual.theta << ',' << actual.rho << '\n';
    } else if (!close(actual.vega, row.expected.vega) || !close(actual.theta, row.expected.theta)
               || !close(actual.rho, row.expected.rho)) {
        std::cerr << std::setprecision(17) << method << ": got " << actual.vega << ", "
                  << actual.theta << ", " << actual.rho << " expected " << row.expected.vega
                  << ", " << row.expected.theta << ", " << row.expected.rho << '\n';
        throw std::runtime_error(row.name + "/" + row.kind);
    }
}

template <class Error, class Function>
void expect_error(Function fn) {
    try {
        fn();
    } catch (const Error&) {
        return;
    }
    throw std::runtime_error("Expected exception was not thrown");
}

void check_contract() {
    const std::vector<double> empty, one{1};
    const CallInputs no_rows{empty, empty, empty, empty, empty, empty};
    const auto result = evaluate_greeks(no_rows, OptionKind::Put, BatchBackend::ParUnseq);
    if (!result.vega.empty() || !result.theta.empty() || !result.rho.empty()) {
        throw std::runtime_error("Empty input returned values");
    }
    expect_error<std::invalid_argument>([] { evaluate_greeks(-1, 1, 0, 0, 1, 1); });
    expect_error<std::invalid_argument>([] { evaluate_greeks(1, 1, 0, 0, 1, -1); });
    expect_error<std::invalid_argument>([] { evaluate_greeks(1, 1, 0, 0, 1, 1, static_cast<OptionKind>(42)); });
    expect_error<std::invalid_argument>([&] {
        evaluate_greeks(no_rows, OptionKind::Call, static_cast<BatchBackend>(42));
    });
    expect_error<std::invalid_argument>([&] { evaluate_greeks(CallInputs{one, empty, one, one, one, one}); });
    expect_error<std::overflow_error>([] { evaluate_greeks(100, 100, -1000, -1000, 0.2, 1); });
    const std::vector<double> rate{-1000};
    expect_error<std::overflow_error>([&] {
        evaluate_greeks(CallInputs{one, one, rate, rate, one, one}, OptionKind::Call, BatchBackend::ParUnseq);
    });
}
} // namespace

int main(int argc, char** argv) {
    try {
        const bool dump = argc == 3 && std::string(argv[2]) == "--dump";
        if (argc != 2 && !dump) throw std::runtime_error("Usage: bs_greeks_contract_test CASES.csv [--dump]");
        const auto rows = read_cases(argv[1]);
        if (dump) std::cout << std::setprecision(17) << "case,kind,method,vega,theta,rho\n";
        for (const auto kind : {OptionKind::Call, OptionKind::Put}) {
            std::vector<const Row*> selected;
            std::array<std::vector<double>, 6> inputs;
            for (const auto& row : rows) {
                if ((row.kind == "call") != (kind == OptionKind::Call)) continue;
                const auto& a = row.inputs;
                check(row, "cpp_scalar", evaluate_greeks(a[0], a[1], a[2], a[3], a[4], a[5], kind), dump);
                selected.push_back(&row);
                for (std::size_t j = 0; j < a.size(); ++j) inputs[j].push_back(a[j]);
            }
            const CallInputs batch{inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5]};
            for (const auto backend : {BatchBackend::Serial, BatchBackend::ParUnseq}) {
                const char* method = backend == BatchBackend::Serial ? "cpp_batch" : "cpp_par_unseq";
                const auto actual = evaluate_greeks(batch, kind, backend);
                for (std::size_t i = 0; i < selected.size(); ++i) {
                    check(*selected[i], method, {actual.vega[i], actual.theta[i], actual.rho[i]}, dump);
                }
            }
        }
        if (!dump) {
            check_contract();
            std::cout << rows.size() << " golden rows (scalar, serial and par_unseq) and contract checks passed\n";
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
