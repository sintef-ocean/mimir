#!/usr/bin/env bash
if [ -z ${OSPL_HOME+x} ]; then
    source /opt/opensplice-rts/release.com;
else
    case ":$PATH:" in
        *:${OSPL_HOME}/bin:*);;
        *) export PATH=${OSPL_HOME}/bin${PATH:+:}${PATH};;
    esac
fi
case ":$PATH:" in
    *:/opt/bin:*);;
    *) export PATH=/opt/bin${PATH:+:}${PATH};;
esac
case ":$PATH:" in
  *:/opt/fkin/bin:*);;
  *) export PATH=/opt/fkin/bin${PATH:+:}${PATH};;
esac
