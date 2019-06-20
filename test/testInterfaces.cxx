#include <iostream>
#include <vectorFlow/PipelineFlow.h>
#include <vectorFlow/ComplexFlow.h>
#include <vectorFlow/RngWrapper.h>

using namespace vectorflow;

/// The user state
struct MyState {
  int fId = 0;    ///< unique id
  int fCount = 0; ///< How many processing steps for this state
  double fVal;    ///< value

  MyState(int id, double val) : fId(id), fVal(val) {}
  void Print() { std::cout << "  (" << fId << ":  " << fVal << ", count=" << fCount << ")\n"; }
};

//_______________________________________________________________________________________________
/// The user task A. Dispatches state to multiplier client state < 9, to divider client if state > 10
/// to averaging client if in interval
struct TaskA : public Work<MyState, std::vector<MyState*>>
{  
  static constexpr size_t kMultiplierClient = 0;
  static constexpr size_t kDividerClient    = 1;
  static constexpr size_t kAveragingClient  = 2;

  /// Dispatcher function
  void SelectClient(MyState * state) {
    if (state->fVal < 9)
      Dispatch(state, kMultiplierClient);
    else if (state->fVal > 10)
      Dispatch(state, kDividerClient);
    else
      Dispatch(state, kAveragingClient);    
  }

  /// Scalar executor
  void Execute(MyState * state) {
    std::cout << "executing TaskSelector scalar\n";    
    SelectClient(state);
    state->fCount++;
    state->Print();
  }

  // Vector executor
  void Execute(std::vector<MyState*> const &states)
  {
    std::cout << "executing TaskSelector vector\n";
    for (auto state : states) {
      SelectClient(state);
      state->fCount++;
      state->Print();
    }
  }
};

//_______________________________________________________________________________________________
/// The user task B. Multiplies the state value by 2
struct TaskB : public Work<MyState, std::vector<MyState*>>
{
  static constexpr size_t kSelectoClient = 0;

  /// Scalar executor
  void Execute(MyState * state) {
    std::cout << "executing TaskMultiplier scalar\n";
    state->fVal *= 2;
    state->fCount++;
    state->Print();
    Dispatch(state, kSelectoClient);
  }

  // Vector executor
  void Execute(std::vector<MyState*> const &states)
  {
    // Normally this code should be vectorized on the states vector
    std::cout << "executing TaskMultiplier vector\n";
    for (auto state : states) {
      state->fVal *= 2;
      state->fCount++;
      state->Print();
      Dispatch(state, kSelectoClient);
    }
  }
};

//_______________________________________________________________________________________________
/// The user task C. Divides the state value by 1.5
struct TaskC : public Work<MyState, std::vector<MyState*>>
{
  static constexpr size_t kSelectoClient = 0;

  /// Scalar executor
  void Execute(MyState * state) {
    std::cout << "executing TaskDivider scalar\n";
    state->fVal /= 1.5;
    state->fCount++;
    state->Print();
    Dispatch(state, kSelectoClient);
  }

  // Vector executor
  void Execute(std::vector<MyState*> const &states)
  {
    // Normally this code should be vectorized on the states vector
    std::cout << "executing TaskDivider vector\n";
    for (auto state : states) {
      state->fVal /= 1.5;
      state->fCount++;
      state->Print();
      Dispatch(state, kSelectoClient);
    }
  }
};

//_______________________________________________________________________________________________
/// The user task D. Does the average of the processed state values
struct TaskD : public Work<MyState, std::vector<MyState*>>
{
  size_t fNstates = 0;
  double fSum = 0.;

  double Average() const { return (fNstates) ? fSum / fNstates : 0.; }

  /// Scalar executor
  void Execute(MyState * state) {
    std::cout << "executing TaskAveraging scalar\n";
    fNstates++;
    fSum += state->fVal;
    state->fCount++;
    state->Print();
  }

  // Vector executor
  void Execute(std::vector<MyState*> const &states)
  {
    std::cout << "executing TaskAveraging vector\n";
    for (auto state : states) {
      fNstates++;
      fSum += state->fVal;
      state->fCount++;
      state->Print();
    }
  }
};


//_______________________________________________________________________________________________
/// Basic test for interfaces
int main(int argc, char const *argv[])
{
  using MyStateFlow_t = ComplexFlow<MyState, std::vector<MyState*>, 4>;
  
  static constexpr size_t kSelectorStage   = 0;
  static constexpr size_t kMultiplierStage = 1;
  static constexpr size_t kDividerStage    = 2;
  static constexpr size_t kAveragingStage  = 3;

  TaskA tA;
  TaskB tB;
  TaskC tC;
  TaskD tD;

  /// Make the following flow:
  ///
  ///          taskB
  ///        //     
  ///   taskA ----- taskD  
  ///        \\     
  ///          taskC
  ///
  /// TaskA takes a state (number). If the number is less than 9 it sends the state to taskB, if
  /// greater than 10 to taskC, otherwise to taskD. Both task B and C alter the state, then send 
  /// it back to taskA for selection.

  /// Add all tasks in the flow
  MyStateFlow_t convergenceFlow;
  convergenceFlow.AddWork(&tA, kSelectorStage);
  convergenceFlow.AddWork(&tB, kMultiplierStage);
  convergenceFlow.AddWork(&tC, kDividerStage);
  convergenceFlow.AddWork(&tD, kAveragingStage);

  // We want taskB and taskC to run in vector mode
  convergenceFlow.SetVectorMode(kMultiplierStage, true);
  convergenceFlow.SetVectorMode(kDividerStage, true);

  // Get references to the input baskets of the different tasks
  std::vector<MyState*> &bSelector   = convergenceFlow.InputData(kSelectorStage);
  std::vector<MyState*> &bMultiplier = convergenceFlow.InputData(kMultiplierStage);
  std::vector<MyState*> &bDivider    = convergenceFlow.InputData(kDividerStage);
  std::vector<MyState*> &bAveraging  = convergenceFlow.InputData(kAveragingStage);
  
  /// Connect the clients to the tasks (see diagram above)
  tA.AddClient(&bMultiplier);
  tA.AddClient(&bDivider);
  tB.AddClient(&bSelector);
  tC.AddClient(&bSelector);

  /// Produce some random data between [0, 100]
  vectorflow::RngWrapper rng;
  std::vector<MyState> input;
  for (auto i = 0; i < 100; ++i) {
    input.push_back(MyState(i, rng.uniform(0., 100.)));
    convergenceFlow.AddData(&input[i]); // we need to add pointers to existing objects
  }

  /// Process the complex flow as long as there is still data in the buffers
  while (convergenceFlow.GetNstates())
    convergenceFlow.Execute();

  return 0;
}
