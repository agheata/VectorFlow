#ifndef VECTORFLOW_WORK_H
#define VECTORFLOW_WORK_H

namespace vectorflow {

/// Base class for generic work aware of the execution mode: scalar or vector.
/** Must be constructed with a scalar data type and a container type for pointers of this data, e.g:
      Work<Track, std::vector<Track*>>
  The work may optionally dispatch the input state to a follow0up work task, via a DataContainer. To do this,
  one should first register the input containers of the client tasks, then call the Dispatch method to the 
  appropriate one from the implementation of the Execute methods.
*/
template <typename Data, typename DataContainer>
struct Work
{
  using DataPtr_t     = Data*;
  using ClientsList_t = std::vector<DataContainer*>;

  ClientsList_t fClients;

  virtual ~Work() {};

  /// Add a container of a client work task
  inline void  AddClient(DataContainer *client) { fClients.push_back(client); }

  /// Dispatch current state to a client work task input
  inline void  Dispatch(DataPtr_t state, size_t client) { fClients[client]->push_back(state); }

  /// Dispatch all the input states to a client
  inline void  DispatchAll(DataContainer const &input, size_t client)
  {
    std::copy(input.begin(), input.end(), std::back_inserter(*fClients[client]));
  }

  virtual void Execute(DataPtr_t)             = 0;
  virtual void Execute(DataContainer const &) = 0;

  inline void ExecuteLoop(DataContainer const &input) { for ( DataPtr_t dataptr : input) Execute(dataptr); }
};
} // namespace vectorflow
#endif
