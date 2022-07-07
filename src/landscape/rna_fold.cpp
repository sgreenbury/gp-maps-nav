#include "rna_fold.hpp"

std::string GetProjectHome(const char *key = "GP_MAPS_NAV") {
  char *project_home = std::getenv(key);
  if (project_home == nullptr) {
    throw std::runtime_error("Project root path not defined");
  }
  return std::string(project_home);
}

void Fold(char *&char_seq, char *&phys_seq) { fold(char_seq, phys_seq); }

float InverseFold(char *&char_seq, char *&phys_seq) {
  return inverse_fold(char_seq, phys_seq);
}

// Function taken from: https://stackoverflow.com/a/478960
std::string exec(const char *cmd) {
  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  return result;
}

std::string shape(std::string in, int level = 0) {
  if (level == 0) {
    return in;
  }
  std::string rnashapes_path = "RNAshapes";
  std::string run_str =
      absl::StrFormat("%s -D\"%s\" -t %d", rnashapes_path, in, level);
  std::string shape_str = exec(run_str.c_str());
  if (shape_str.back() == '\n') {
    return shape_str.erase(shape_str.size() - 1);
  }
  return shape_str;
}
