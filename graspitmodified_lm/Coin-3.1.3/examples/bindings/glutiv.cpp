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

// glutiv.cpp -- example demonstrating Coin (or Open Inventor) bound
// to the GLUT user interface abstraction.
//

// If you have Coin and the GLUT library properly installed, you should
// be able to build by simply doing:
//
// 	$ coin-config --build glutiv glutiv.cpp -lglut
//
// (or -lglut32 if you're on MSWindows with Cygwin).
//
// Note that to compile on Mac OS X, you have to link against the
// GLUT framework, like so:
// 
//  $ export LDFLAGS="-framework GLUT -lobjc"
//  $ coin-config --build-app glutiv glutiv.cpp 


// *************************************************************************

// FIXME: note that there are some limitations to this example:
//
//  * The most important one is that events are not translated from
//    native X11 events to Coin events (to be sent to the scene
//    graph). This means that for instance draggers and manipulators in
//    the scene graph will not respond to attempts at interaction.
//
//  * The sensor queue processing is just a hack. Don't use this
//    in production code.
//
// 20031113 mortene.

// *************************************************************************

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif 

#include <Inventor/SoDB.h>
#include <Inventor/SoInteraction.h>
#include <Inventor/SoSceneManager.h>
#include <Inventor/nodekits/SoNodeKit.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoRotor.h>
#include <Inventor/nodes/SoShuttle.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoTransform.h>

// *************************************************************************

const int SCENEWINDOWS = 3;

SoSceneManager * scenemanager[SCENEWINDOWS];
int glutwin[SCENEWINDOWS];

// *************************************************************************

int
winid2idx(const int winid)
{
  for (int i=0; i < SCENEWINDOWS; i++) {
    if (winid == glutwin[i]) return i;
  }
}

// Redraw on scenegraph changes.
void
redraw_cb(void * user, SoSceneManager * manager)
{
  unsigned int idx = (uintptr_t)user;

  glutSetWindow(glutwin[idx]);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  scenemanager[idx]->render();

  glutSwapBuffers();
}

// Redraw on expose events.
void
expose_cb(void)
{
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  scenemanager[winid2idx(glutGetWindow())]->render();

  glutSwapBuffers();
}

// Reconfigure on changes to window dimensions.
void
reshape_cb(int w, int h)
{
  int idx = winid2idx(glutGetWindow());

  scenemanager[idx]->setWindowSize(SbVec2s(w, h));
  scenemanager[idx]->setSize(SbVec2s(w, h));
  scenemanager[idx]->setViewportRegion(SbViewportRegion(w, h));
  scenemanager[idx]->scheduleRedraw();
}

// Process the internal Coin queues when idle. Necessary to get the
// animation to work.
void
idle_cb(void)
{
  SoDB::getSensorManager()->processTimerQueue();
  SoDB::getSensorManager()->processDelayQueue(TRUE);
}

// *************************************************************************

// Make the common scenegraph. Just a little silly something to show
// off some geometry nodes and simple animation.

SoSeparator *
commongraph(void)
{
  SoSeparator * root = NULL;

  if (!root) {
    root = new SoSeparator;

    {
      SoSeparator * cylsep = new SoSeparator;
      root->addChild(cylsep);

      SoRotor * rotor = new SoRotor;
      rotor->speed.setValue(0.1f);
      cylsep->addChild(rotor);

      cylsep->addChild(new SoCylinder);
    }

    SoMaterial * material = new SoMaterial;
    material->diffuseColor.setValue(0.0f, 0.0f, 0.0f);
    root->addChild(material);

    SoSphere * sphere = new SoSphere;

    SoTransform * trans0 = new SoTransform;
    trans0->translation.setValue(-0.5f, 0.5f, 1.0f);
    trans0->scaleFactor.setValue(0.1f, 0.1f, 0.1f);
    root->addChild(trans0);

    root->addChild(sphere);

    SoShuttle * shuttle = new SoShuttle;
    shuttle->translation0.setValue(0.0f, 5.0f, 0.0f);
    shuttle->translation1.setValue(0.0f, -5.0f, 0.0f);
    shuttle->on.setValue(TRUE);
    root->addChild(shuttle);

    SoTransform * trans1 = new SoTransform;
    trans1->translation.setValue(10.0f, 0.0f, 0.0f);
    root->addChild(trans1);

    root->addChild(sphere);
  }

  return root;
}

// *************************************************************************

#ifdef _WIN32

#include <windows.h>
#include <winbase.h>

int
WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine,
        int nCmdShow)

#else // UNIX

int
main(int argc, char ** argv)

#endif

{
  // initialize Coin and glut libraries

  SoDB::init();
  SoNodeKit::init();
  SoInteraction::init();

#ifdef _WIN32
  int argc = 1;
  char * argv[] = { "glutiv.exe", (char *) NULL };
  glutInit(&argc, argv);
#else
  glutInit(&argc, argv);
#endif
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

  // Note: _don't_ use SoGroup, as TGS' Inventor has a bug where
  // lightsource contribution will get accumulated over runs.
  SoSeparator * root[SCENEWINDOWS];

  for (unsigned int i=0; i < SCENEWINDOWS; i++) {
    // set up individual parts of scenegraph

    root[i] = new SoSeparator;
    root[i]->ref();
    SoPerspectiveCamera * camera = new SoPerspectiveCamera;
    root[i]->addChild(camera);
    root[i]->addChild(new SoDirectionalLight);
    SoDrawStyle * drawstyle = new SoDrawStyle;
    drawstyle->style.setValue(i % 3);
    root[i]->addChild(drawstyle);
    root[i]->addChild(commongraph());

    // initialize scenemanager instance

    scenemanager[i] = new SoSceneManager;
    scenemanager[i]->setRenderCallback(redraw_cb, (void *)i);
    scenemanager[i]->setBackgroundColor(SbColor(0.2f, 0.2f, i * 0.2f));
    scenemanager[i]->activate();
    camera->viewAll(root[i], scenemanager[i]->getViewportRegion());
    scenemanager[i]->setSceneGraph(root[i]);

    // glut window initialization

    glutInitWindowSize(512, 400);
    SbString title("window ");
    title += (char)(i + 0x30);
    glutwin[i] = glutCreateWindow(title.getString());
    glutDisplayFunc(expose_cb);
    glutReshapeFunc(reshape_cb);
  }

  SoDB::setRealTimeInterval(1/120.0);

  // start main loop processing (with an idle callback)

  glutIdleFunc(idle_cb);
  glutMainLoop();

  // clean up Coin resource use

  for (int j=0; j < SCENEWINDOWS; j++) {
    root[j]->unref();
    delete scenemanager[j];
  }

  return 0;
}

// *************************************************************************
