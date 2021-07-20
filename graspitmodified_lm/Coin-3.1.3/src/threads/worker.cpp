/**************************************************************************\
 *
 *  This file is part of the Coin 3D visualization library.
 *  Copyright (C) by Kongsberg Oil & Gas Technologies.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  ("GPL") version 2 as published by the Free Software Foundation.
 *  See the file LICENSE.GPL at the root directory of this source
 *  distribution for additional information about the GNU GPL.
 *
 *  For using Coin with software that can not be combined with the GNU
 *  GPL, and for taking advantage of the additional benefits of our
 *  support services, please contact Kongsberg Oil & Gas Technologies
 *  about acquiring a Coin Professional Edition License.
 *
 *  See http://www.coin3d.org/ for more information.
 *
 *  Kongsberg Oil & Gas Technologies, Bygdoy Alle 5, 0257 Oslo, NORWAY.
 *  http://www.sim.no/  sales@sim.no  coin-support@coin3d.org
 *
\**************************************************************************/

#include <Inventor/C/threads/worker.h>

#include <stdlib.h>
#include <assert.h>

#include <Inventor/C/threads/thread.h>
#include <Inventor/C/threads/mutex.h>
#include <Inventor/C/threads/condvar.h>
#include <Inventor/C/errors/debugerror.h>

#include "threads/workerp.h"
#include "threads/mutexp.h"

/* ********************************************************************** */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
 * this is the thread's main loop.
 */
static void 
worker_thread_loop(cc_worker * worker)
{
  cc_mutex_lock(worker->mutex);
  cc_mutex_lock(worker->beginmutex);
  /* signal master that we've entered the scheduling loop */
  cc_condvar_wake_one(worker->begincond);
  cc_mutex_unlock(worker->beginmutex);

  while (!worker->shutdown) {
    /* wait for job */
    cc_condvar_wait(worker->cond, worker->mutex);
    /* check if there are more jobs */
    if (!worker->shutdown) {
      /* do the scheduled job */
      worker->workfunc(worker->workclosure);
      /* call the idle cb */
      if (worker->idlecb) {
        worker->idlecb(worker, worker->idleclosure);
      }
    }
  }

  worker->workfunc = NULL;
  /* remember to unlock mutex after we break out of the loop */
  cc_mutex_unlock(worker->mutex);
}

/*
 * entry point for thread. It just jumps to the main loop
 */
static void * 
worker_thread_entry(void * data)
{
  cc_worker * worker = (cc_worker*) data;
  worker_thread_loop(worker);
  return NULL;
}

/*
 * called to start thread. Assumes worker->mutex is locked by caller */
static void 
worker_start_thread(cc_worker * worker)
{
  if (!worker->threadisrunning) {
    cc_mutex_lock(worker->beginmutex);
    /* Unlock worker mutex before starting thread. The new thread will use
       mutex and beginmutex  to synchronize */
    cc_mutex_unlock(worker->mutex);
    worker->thread = cc_thread_construct(worker_thread_entry, worker);
    
    /* Wait for thread to get to the main loop. The new thread will
     have worker->mutex locked when this signal arrives */
    cc_condvar_wait(worker->begincond, worker->beginmutex);
    
    /* lock thread again so that we know the new thread is waiting for
       a signal before we return */
    cc_mutex_lock(worker->mutex);
    worker->threadisrunning = 1;
    cc_mutex_unlock(worker->beginmutex);
  }
}

/*
 * join caller thread with this thread.
 */
static void 
worker_stop_thread(cc_worker * worker)
{
  if (worker->threadisrunning) {
    cc_mutex_lock(worker->mutex);
    worker->threadisrunning = 0;
    
    worker->shutdown = TRUE;  /* signal thread to exit loop */
    /* in case thread is waiting for cond... */
    cc_condvar_wake_one(worker->cond);
    cc_mutex_unlock(worker->mutex);
    /* wait for thread to finish */
    cc_thread_join(worker->thread, NULL);
    cc_thread_destruct(worker->thread);
    worker->thread = NULL;
    worker->shutdown = FALSE; /* reset signal */
  }
}
  

cc_worker * 
cc_worker_construct(void)
{
  cc_worker * worker = (cc_worker*) malloc(sizeof(cc_worker));
  assert(worker);

  worker->mutex = cc_mutex_construct();
  worker->cond = cc_condvar_construct();
  worker->begincond = cc_condvar_construct();
  worker->beginmutex = cc_mutex_construct();
  worker->thread = NULL; /* delay creating thread */
  worker->threadisrunning = FALSE;
  worker->shutdown = FALSE;
  worker->workfunc = NULL;
  worker->workclosure = NULL;
  worker->idlecb = NULL;
  worker->idleclosure = NULL;
  return worker;
}

/*!
  Will wait for the current task to be finished before destructing the worker.
*/
void 
cc_worker_destruct(cc_worker * worker)
{
  if (worker->threadisrunning) {
    worker_stop_thread(worker);
  }
  assert(worker->thread == NULL);
  cc_mutex_destruct(worker->mutex);
  cc_condvar_destruct(worker->cond);
  cc_condvar_destruct(worker->begincond);
  cc_mutex_destruct(worker->beginmutex);
  free(worker);
}

/*!
  Start the worker by either waking is from a condvar_wait() or
  creating a new thread if the thread is not running.
*/
SbBool 
cc_worker_start(cc_worker * worker, cc_worker_f * workfunc, void * closure)
{
  assert(workfunc);

  cc_mutex_lock(worker->mutex);  
  worker->workfunc = workfunc;
  worker->workclosure = closure;

  if (!worker->threadisrunning) {
    worker_start_thread(worker);
  }
  
  /* We now know that thread is waiting for a signal */
  cc_condvar_wake_one(worker->cond);
  cc_mutex_unlock(worker->mutex);
  return TRUE;
}

SbBool 
cc_worker_is_busy(cc_worker * worker)
{
  SbBool busy = TRUE;
  
  if (cc_mutex_try_lock(worker->mutex)) {
    busy = FALSE;
    cc_mutex_unlock(worker->mutex);
  }
  return busy;
}

void 
cc_worker_wait(cc_worker * worker)
{
  /* synchronize by stopping thread. A new thread must be created
   * if more work needs to be done, but this is ok, I guess.
   */
  worker_stop_thread(worker);
}

void 
cc_worker_set_idle_callback(cc_worker * worker, cc_worker_idle_f * cb, 
                            void * closure)
{
  cc_mutex_lock(worker->mutex);
  worker->idlecb = cb;
  worker->idleclosure = closure;
  cc_mutex_unlock(worker->mutex);
}

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */
