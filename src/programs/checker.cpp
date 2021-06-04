#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>

#include <yaml-cpp/yaml.h>
#include <casadi/casadi.hpp>

#include <vector>
#include <string>
#include <iostream>
//#include <tuple>

#include <mimir/program/Config.hpp>


int main()
{

  // To be read from a config file
  YAML::Node config = YAML::LoadFile("./integrator.yml");
  YAML::Node schema = YAML::LoadFile("./integrator_schema.yml");

  YAML::Node combo;
  combo["config"] = config;
  combo["schema"] = schema;

  try {
    auto dict = mimir::program::parse_config(combo["config"], combo["schema"]);

    using std::vector;
    using namespace casadi;

     // Ode system
    MX x = MX::sym("x", 4);
    MX A = MX(4, 4);
    A(0,1) = 1;
    A(1,1) = -.1;
    A(2,3) = 1;
    A(3,3) = -.1;
    MX rhs = mtimes(A,x);

    MXDict ode;
    ode["x"] = x;
    ode["ode"] = rhs;

    // Time grid output
    auto steps = static_cast<size_t>(dict["t_final"].as_double()/dict["t_step"].as_double());
    vector<double> grid;
    for (size_t i = 0; i < steps; ++i)
      grid.push_back(i*dict["t_step"].as_double());

    Dict generic;
    generic["grid"] = grid;
    generic["output_t0"] = true;

    Dict grator = dict["integrator"];
    Dict combo2 = combine(grator["settings"], generic, true);

    std::cout << doc_integrator(grator["name"].as_string()) << std::endl;
    // Integrator
    auto F = integrator("F", grator["name"].as_string(), ode, combo2);

    // Evaluate F
    DMDict res = F(DMDict{{"x0", dict["x0"].as_double_vector()}});

    std::cout << "Result with " << grator["name"] << ": " << res["xf"] << std::endl;
  }
  catch(YAML::Exception&e)
  {
    std::cout << e.what() << std::endl;
  }


  return 0;
}
