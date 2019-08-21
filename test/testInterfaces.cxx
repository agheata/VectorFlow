#include <iostream>
#include <vectorFlow/vectorFlow>

using namespace vectorflow;

bool kVerbose = false;

/// The user state
struct MyState {
  int fId = 0;    ///< unique id
  int fCount = 0; ///< How many processing steps for this state
  double fVal;    ///< value

  MyState(int id, double val) : fId(id), fVal(val) {}
  void Print() {
    std::cout << "  (" << fId << ":  " << fVal << ", count=" << fCount << ")\n";
  }
};

//_______________________________________________________________________________________________
/// The user task A. Dispatches state to multiplier client state < 9, to divider
/// client if state > 10 to averaging client if in interval
struct TaskSelector : public Work<MyState, std::vector<MyState *>> {
  static constexpr size_t kMultiplierClient = 0;
  static constexpr size_t kDividerClient = 1;
  static constexpr size_t kAveragingClient = 2;

  /// Dispatcher function
  void SelectClient(MyState *state) {
    if (state->fVal < 9)
      Dispatch(state, kMultiplierClient);
    else if (state->fVal > 10)
      Dispatch(state, kDividerClient);
    else
      Dispatch(state, kAveragingClient);
  }

  /// Scalar executor
  void Execute(MyState *state) {
    SelectClient(state);
    state->fCount++;
  }

  // Vector executor
  void Execute(std::vector<MyState *> const &states) {
    if (kVerbose)
      std::cout << "TaskSelector::Execute not implemented for vector mode\n";
  }
};

//_______________________________________________________________________________________________
/// The user task B. Multiplies the state value by 2
struct TaskMultiplier : public Work<MyState, std::vector<MyState *>> {
  static constexpr size_t kSelectorClient = 0;
  size_t fNcalls = 0;
  size_t fVecSize = 0;

  /// Scalar executor
  void Execute(MyState *state) {
    state->fVal *= 2;
    state->fCount++;
    Dispatch(state, kSelectorClient);
  }

  // Vector executor
  void Execute(std::vector<MyState *> const &states) {
    // Normally this code should be vectorized on the states vector
    if (kVerbose)
      std::cout << "executing TaskMultiplier vector on " << states.size()
                << " states\n";
    fNcalls++;
    fVecSize += states.size();
    for (auto state : states) {
      state->fVal *= 2;
      state->fCount++;
      Dispatch(state, kSelectorClient);
    }
  }
};

//_______________________________________________________________________________________________
/// The user task C. Divides the state value by 1.5
struct TaskAveragingivider : public Work<MyState, std::vector<MyState *>> {
  static constexpr size_t kSelectorClient = 0;
  size_t fNcalls = 0;
  size_t fVecSize = 0;

  /// Scalar executor
  void Execute(MyState *state) {
    state->fVal /= 1.5;
    state->fCount++;
    Dispatch(state, kSelectorClient);
  }

  // Vector executor
  void Execute(std::vector<MyState *> const &states) {
    // Normally this code should be vectorized on the states vector
    if (kVerbose)
      std::cout << "executing TaskAveragingivider vector on " << states.size()
                << " states\n";
    fNcalls++;
    fVecSize += states.size();
    for (auto state : states) {
      state->fVal /= 1.5;
      state->fCount++;
      Dispatch(state, kSelectorClient);
    }
  }
};

//_______________________________________________________________________________________________
/// The user task D. Does the average of the processed state values
struct TaskAveraging : public Work<MyState, std::vector<MyState *>> {
  size_t fNstates = 0;
  double fSum = 0.;

  double Average() const { return (fNstates) ? fSum / fNstates : 0.; }

  /// Scalar executor
  void Execute(MyState *state) {
    fNstates++;
    fSum += state->fVal;
    state->fCount++;
    if (kVerbose) {
      std::cout << "=== final value for state: ";
      state->Print();
    }
  }

  // Vector executor
  void Execute(std::vector<MyState *> const &states) {
    if (kVerbose)
      std::cout
          << "TaskSelectorveraging::Execute not implemented for vector mode\n";
  }
};

//_______________________________________________________________________________________________
/// Basic test for interfaces
int main(int argc, char const *argv[]) {
  using MyStateFlow_t = ComplexFlow<MyState, std::vector<MyState *>, 4>;

  static constexpr size_t kSelectorStage = 0;
  static constexpr size_t kMultiplierStage = 1;
  static constexpr size_t kDividerStage = 2;
  static constexpr size_t kAveragingStage = 3;

  int ninput = 0;
  if (argc == 1) {
    std::cout << "=== Usage: testInterfaces N [v] \n===      will run on N "
                 "states, optionally verbose.\n";
    return 1;
  }
  if (argc > 1)
    ninput = std::atoi(argv[1]);
  if (argc > 2)
    kVerbose = true;

  TaskSelector tA;
  TaskMultiplier tB;
  TaskAveragingivider tC;
  TaskAveraging tD;

  /// Make the following flow:
  ///
  ///          TaskMultiplier
  ///        //
  ///   TaskSelector ----- TaskAveraging
  ///        \\     
  ///          TaskAveragingivider
  ///
  /// TaskSelector takes a state (number). If the number is less than 9 it sends
  /// the state to TaskMultiplier, if greater than 10 to TaskAveragingivider,
  /// otherwise to TaskAveraging. Both task B and C alter the state, then send
  /// it back to TaskSelector for selection. All values of states reaching
  /// TaskAveraging have to fall in [9, 10] and their computed average as well.
  /// The goal is to maximize the input state population for TaskMultiplier and
  /// TaskAveragingivider, which are the tasks providing a vector implementation
  /// of Execute.

  /// Add all tasks in the flow
  MyStateFlow_t convergenceFlow;
  convergenceFlow.AddWork(&tA, kSelectorStage);
  convergenceFlow.AddWork(&tB, kMultiplierStage);
  convergenceFlow.AddWork(&tC, kDividerStage);
  convergenceFlow.AddWork(&tD, kAveragingStage);

  // We want TaskMultiplier and TaskAveragingivider to run in vector mode
  convergenceFlow.SetVectorMode(kMultiplierStage, true);
  convergenceFlow.SetVectorMode(kDividerStage, true);

  // Get references to the input baskets of the different tasks
  std::vector<MyState *> &bSelector = convergenceFlow.InputData(kSelectorStage);
  std::vector<MyState *> &bMultiplier =
      convergenceFlow.InputData(kMultiplierStage);
  std::vector<MyState *> &bDivider = convergenceFlow.InputData(kDividerStage);
  std::vector<MyState *> &bAveraging =
      convergenceFlow.InputData(kAveragingStage);

  /// Connect the clients to the tasks (see diagram above)
  tA.AddClient(&bMultiplier);
  tA.AddClient(&bDivider);
  tA.AddClient(&bAveraging);
  tB.AddClient(&bSelector);
  tC.AddClient(&bSelector);

  /// Produce some random data between [0, 100]
  vectorflow::RngWrapper rng;
  std::vector<MyState> input;
  std::cout << "Executing testInterfaces on " << ninput << " states\n";
  for (auto i = 0; i < ninput; ++i)
    input.push_back(MyState(i, rng.uniform(0., 100.)));

  /// Add the inpit data to the flow
  for (auto &state : input)
    convergenceFlow.AddData(&state);

  /// Process the complex flow as long as there is still data in the buffers
  while (convergenceFlow.GetNstates())
    convergenceFlow.Execute();

  std::cout << "===== Average of values for final states: "
            << tD.fSum / tD.fNstates << std::endl;
  std::cout << "===== TaskMultiplier::Execute had an average input of: "
            << double(tB.fVecSize) / tB.fNcalls << " states\n";

  std::cout << "===== TaskAveragingivider::Execute had an average input of: "
            << double(tC.fVecSize) / tC.fNcalls << " states\n";

  return 0;
}
